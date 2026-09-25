using EngineLayer.GlycoSearch;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace EngineLayer
{
    // in our database, the N-glycan.gdb should be correct to the new format
    // the class for loading glycan database then creeat the glycan object.
    public static class GlycanDatabase 
    {
        
        /// <summary>
        /// Load Glycan from the database file. Generally, glycan-ions should be generated for N-Glycopepitdes which produce Y-ions; MS method couldn't produce o-glycan-ions
        /// </summary>
        /// <param name="filePath"> Database file path</param>
        /// <param name="ToGenerateIons"> Do we need to generate the glycanIon? </param>
        /// <param name="IsOGlycanSearch"></param>
        /// <returns> A glycan object collection </returns>
        public static IEnumerable<Glycan> LoadGlycan(string filePath, bool ToGenerateIons, bool IsOGlycan)
        {
            // The format is inferred from the first DATA line; it is never declared. Comment and blank lines
            // are skipped while sniffing because a documented database -- such as the seeded custom one --
            // opens with a '#' banner, and judging the format from that banner routes the whole file to the
            // wrong parser. A file with no data line is a legal empty database: either parser yields nothing
            // for it, so which one is chosen does not matter.
            bool isKind = true;
            using (StreamReader lines = new StreamReader(filePath))
            {
                while(lines.Peek() != -1)
                {
                    string line = lines.ReadLine();
                    if (IsCommentOrBlank(line))
                    {
                        continue;
                    }
                    // Structure lines are nested parentheses and always open with one; composition lines
                    // are name-and-count and never do. This is the same test ValidateGlycanLine and
                    // FormatOfExistingEntries use, and it has to be the same one: the validator tells the
                    // user their entry was accepted, so anything it calls a composition has to load as a
                    // composition. The older test -- "does the line contain HexNAc" -- disagreed for every
                    // composition without that literal substring, so a validated Hex(1) was written and
                    // then read back as a structure, where 'e' is not a monosaccharide code. Classification
                    // is unchanged for every database that ships and every test fixture: the structure ones
                    // open with '(', the composition ones open with "HexNAc(".
                    if (line.TrimStart().StartsWith("(", StringComparison.Ordinal))
                    {
                        isKind = false;
                    }
                    break;
                }
            }

            if (isKind)
            {
                return LoadKindGlycan(filePath, ToGenerateIons, IsOGlycan); // open the file of the kind format, example: HexNAc(2)Hex(5)NeuAc(1)Fuc(1)
            }
            else
            {
                return LoadStructureGlycan(filePath, IsOGlycan);            // open the file of the structure format, example: (N(H(A))(A))
            }
        }

        /// <summary>
        /// Every glycan in the given databases, without ions, reporting a database that cannot be read
        /// through <paramref name="warn"/> and carrying on with the rest instead of throwing.
        /// </summary>
        /// <remarks>
        /// For the startup load in <see cref="GlobalVariables.SetUpGlobalVariables"/>, which is only
        /// there to give MetaDraw the glycan names. A throw from there stops MetaMorpheus opening, and a
        /// user's own database is hand-edited, so one typo would lock them out of the window that could
        /// fix it. A database that fails here is left in the list the GlycoSearch task offers: a search
        /// that selects it reads it again and stops with the same named error, which is where it matters.
        /// </remarks>
        public static List<Glycan> LoadGlycansOrWarn(IEnumerable<string> databasePaths, bool isOGlycan, Action<string> warn)
        {
            var glycans = new List<Glycan>();
            foreach (string path in databasePaths)
            {
                try
                {
                    // Materialised here: the loaders are lazy, so the parse error surfaces on enumeration,
                    // and a half-read database is not added.
                    glycans.AddRange(LoadGlycan(path, false, isOGlycan).ToList());
                }
                // Anything, not only the loaders' own named errors: a file locked by an editor, or a parse
                // failure deeper in Struct2Glycan, would stop startup just the same.
                catch (Exception ex)
                {
                    warn($"The glycan database '{Path.GetFileName(path)}' could not be read, so its glycans are not " +
                        $"available to MetaDraw in this session. {ex.Message} Fix the file at {path} and restart MetaMorpheus.");
                }
            }
            return glycans;
        }

        public const string MonoSaccharidesHeader = "Name\tSingleCharCode\tMonoisotopicMass\tDiagnosticIonMasses\tDescription";

        /// <summary>
        /// Load custom monosaccharide definitions from a tab-separated file and register them with
        /// Glycan so they are recognized by both glycan-database formats and the structure validator.
        ///
        /// File format (tab-separated, one entry per line, header optional):
        ///   Name  SingleCharCode  MonoisotopicMass  DiagnosticIonMasses  Description
        ///
        /// Lines starting with '#' and blank lines are skipped. The header row (a line beginning
        /// with "Name" followed by a tab) is also skipped if present. Description column is
        /// optional and ignored. DiagnosticIonMasses is optional and may contain zero or more
        /// comma-separated decimal m/z values.
        ///
        /// Mass values are decimal Daltons (e.g. "176.03209") and are internally scaled by 1e5 to
        /// match the integer-mass representation used throughout the glycan code. Diagnostic ion
        /// m/z values are stored as supplied (no hydrogen-mass offset).
        ///
        /// A malformed line throws MetaMorpheusException with the file name, line number, raw line
        /// content, and the specific problem (bad char, non-numeric mass).
        ///
        /// A diagnostic ion that collides with an already-claimed oxonium ion is the one problem that
        /// does NOT throw: the ion is dropped, the monosaccharide is loaded without it, and a warning
        /// goes to GlobalVariables.ErrorsReadingMods. DiagnosticIonMasses shipped in 1.1.8 without
        /// this check, and EnsureCustomMonosaccharideFileExists only writes the template when the file
        /// is missing, so an upgrading user's file can already contain one. Throwing would be fatal
        /// rather than corrective -- LoadGlycans runs from SetUpGlobalVariables before the GUI's
        /// InitializeComponent, so the window never opens and "Open mods/data folder" is unreachable.
        /// </summary>
        public static void LoadCustomMonosaccharides(string filePath)
        {
            if (!File.Exists(filePath))
            {
                return; // file is optional
            }

            int lineNumber = 0;
            using (var reader = new StreamReader(filePath))
            {
                while (reader.Peek() != -1)
                {
                    string line = reader.ReadLine();
                    lineNumber++;
                    if (IsCommentOrBlank(line))
                    {
                        continue;
                    }
                    // Skip the optional column-header row.
                    if (line.TrimStart().StartsWith("Name\t", StringComparison.OrdinalIgnoreCase))
                    {
                        continue;
                    }

                    string[] cols = line.Split('\t');
                    if (cols.Length < 3)
                    {
                        throw new MetaMorpheusException(
                            $"Could not parse custom monosaccharide in '{Path.GetFileName(filePath)}' at line {lineNumber}: \"{line}\". Expected at least 3 tab-separated columns (Name, SingleCharCode, MonoisotopicMass).");
                    }

                    string name = cols[0].Trim();
                    string codeStr = cols[1].Trim();
                    string massStr = cols[2].Trim();
                    string ionsStr = cols.Length > 3 ? cols[3].Trim() : string.Empty;

                    if (codeStr.Length != 1)
                    {
                        throw new MetaMorpheusException(
                            $"Could not parse custom monosaccharide in '{Path.GetFileName(filePath)}' at line {lineNumber}: SingleCharCode must be exactly one character, got \"{codeStr}\".");
                    }
                    char code = codeStr[0];

                    if (!double.TryParse(massStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double massDa))
                    {
                        throw new MetaMorpheusException(
                            $"Could not parse custom monosaccharide in '{Path.GetFileName(filePath)}' at line {lineNumber}: MonoisotopicMass \"{massStr}\" is not a valid decimal number.");
                    }
                    int massScaled = (int)Math.Round(massDa * 1E5);

                    int[] ionsScaled = null;
                    if (!string.IsNullOrWhiteSpace(ionsStr))
                    {
                        var parts = ionsStr.Split(new[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                        ionsScaled = new int[parts.Length];
                        for (int i = 0; i < parts.Length; i++)
                        {
                            if (!double.TryParse(parts[i].Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out double ionDa))
                            {
                                throw new MetaMorpheusException(
                                    $"Could not parse custom monosaccharide in '{Path.GetFileName(filePath)}' at line {lineNumber}: DiagnosticIonMasses entry \"{parts[i]}\" is not a valid decimal number.");
                            }
                            ionsScaled[i] = (int)Math.Round(ionDa * 1E5);
                        }
                    }

                    // Screen the ions before registering rather than letting RegisterCustomMonosaccharide
                    // throw: see the note on this method about why a collision must not be fatal here.
                    if (ionsScaled != null)
                    {
                        ionsScaled = KeepUnclaimedDiagnosticIons(ionsScaled, name, filePath, lineNumber);
                    }

                    try
                    {
                        Glycan.RegisterCustomMonosaccharide(name, code, massScaled, ionsScaled);
                    }
                    catch (ArgumentException ex)
                    {
                        throw new MetaMorpheusException(
                            $"Could not register custom monosaccharide in '{Path.GetFileName(filePath)}' at line {lineNumber}: \"{line}\". {ex.Message}",
                            ex);
                    }

                    WarnThatDiagnosticIonsActAsFilterGates(name, ionsScaled);
                }
            }
        }

        /// <summary>
        /// Returns the subset of <paramref name="ionsScaled"/> that no built-in oxonium ion and no
        /// already-registered custom monosaccharide has claimed, warning about each one dropped.
        /// Returns null when nothing survives, which is what RegisterCustomMonosaccharide expects for
        /// "no diagnostic ions" -- the monosaccharide itself still loads and still works in glycan
        /// databases, it just contributes no ions.
        /// </summary>
        private static int[] KeepUnclaimedDiagnosticIons(int[] ionsScaled, string name, string filePath, int lineNumber)
        {
            List<int> kept = new List<int>(ionsScaled.Length);
            foreach (int ionScaled in ionsScaled)
            {
                string collision = Glycan.DescribeDiagnosticIonCollision(ionScaled);
                if (collision == null)
                {
                    kept.Add(ionScaled);
                    continue;
                }

                AddModWarning(
                    $"Custom monosaccharide '{name}' in '{Path.GetFileName(filePath)}' at line {lineNumber}: " +
                    $"diagnostic ion {(double)ionScaled / 1E5:F5} {collision} " +
                    $"That ion was skipped and '{name}' was loaded without it. Remove or correct it in the file to clear this warning.");
            }

            return kept.Count > 0 ? kept.ToArray() : null;
        }

        /// <summary>
        /// Tells the user, once per startup, that this monosaccharide's diagnostic ions act as strict
        /// accept/reject gates. The other half of the upgrade problem: a file that predates this
        /// feature keeps its column-4 values, OxoniumIonFilt defaults to true, and the banner in the
        /// shipped template explaining all this never reaches them precisely because their file
        /// already exists and is therefore never rewritten.
        /// </summary>
        private static void WarnThatDiagnosticIonsActAsFilterGates(string name, int[] ionsScaled)
        {
            if (ionsScaled == null || ionsScaled.Length == 0)
            {
                return;
            }

            string ions = string.Join(", ", ionsScaled.Select(i => ((double)i / 1E5).ToString("F5", CultureInfo.InvariantCulture)));
            AddModWarning(
                $"Custom monosaccharide '{name}' declares diagnostic ion(s) {ions}. With OxoniumIonFilt enabled (the default) " +
                $"these act as strict gates in O-glycan and N+O-glycan searches: a candidate is rejected when one of these ions " +
                $"is observed but the candidate does not contain '{name}', and when the candidate contains '{name}' but the ion " +
                $"is not observed. Uncheck OxoniumIonFilt to score these ions without filtering on them.");
        }

        /// <summary>
        /// Routes a non-fatal load problem to the notifications area. GlobalVariables.LoadModifications
        /// initialises ErrorsReadingMods one line before LoadGlycans, and the GUI drains it in
        /// PrintErrorsReadingMods once its window is up; the null guard is for callers that load
        /// monosaccharides without SetUpGlobalVariables, i.e. tests.
        /// </summary>
        private static void AddModWarning(string message)
        {
            GlobalVariables.ErrorsReadingMods?.Add(message);
        }

        /// <summary>
        /// Registers a custom monosaccharide in memory (via Glycan.RegisterCustomMonosaccharide) and
        /// persists it to MonosaccharidesCustom.tsv, in the exact format LoadCustomMonosaccharides
        /// reads -- column order, header row, invariant-culture decimal formatting all live here so
        /// the two can't drift apart.
        ///
        /// Return contract: throws MetaMorpheusException if the input is invalid or registration
        /// itself fails -- nothing happened, treat this as a hard failure. Returns null if
        /// registration AND the file write both succeeded. 
        /// Returns a non-null, user-facing warning
        /// if registration succeeded but the file write failed -- the monosaccharide is already
        /// usable for the current session even though it wasn't persisted for the next one.
        /// </summary>
        public static string PersistCustomMonosaccharide(string name, string codeText, string formulaText, string massText, string ionsText, string descriptionText) 
        {
            if (string.IsNullOrEmpty(codeText) || codeText.Length != 1) 
            {
                throw new MetaMorpheusException(
                    $"Could not persist custom monosaccharide: SingleCharCode must be exactly one character, got \"{codeText}\".");
            }
            char code = codeText[0];

            double massDa = !string.IsNullOrEmpty(formulaText) ? Chemistry.ChemicalFormula.ParseFormula(formulaText).MonoisotopicMass
                            : double.Parse(massText, NumberStyles.Float, CultureInfo.InvariantCulture);

            if (massDa <= 0 || massDa > 20000)
            {
                throw new MetaMorpheusException(
                    $"Could not persist custom monosaccharide: MonoisotopicMass must be a positive number below 20000 Da, got {massDa}.");
            }
            int massScaled = (int)Math.Round(massDa * 1E5);

            int[] diagnosticIons = null;
            if (!string.IsNullOrEmpty(ionsText))
            {
                var parsedIons = ionsText.Split(',').Select(p => double.Parse(p.Trim(), NumberStyles.Float, CultureInfo.InvariantCulture)).ToArray();
                foreach (double ionDa in parsedIons)
                {
                    if (ionDa <= 0 || ionDa > 20000)
                    {
                        throw new MetaMorpheusException(
                            $"Could not persist custom monosaccharide: DiagnosticIonMasses entry {ionDa} must be a positive number below 20000 Da.");
                    }
                }
                diagnosticIons = parsedIons.Select(ionDa => (int)Math.Round(ionDa * 1E5)).ToArray();
            }
            try
            {
                Glycan.RegisterCustomMonosaccharide(name, code, massScaled, diagnosticIons);
            }
            catch (ArgumentException ex)
            {
                throw new MetaMorpheusException($"Could not register custom monosaccharide: {ex.Message}", ex);
            }

            //registration already succeeded -- the sugar works this session no matter what happens
            //persist to file so the monosaccharide is still recognized on the next launch
            string line = string.Join("\t", name, code.ToString(), massDa.ToString(CultureInfo.InvariantCulture), ionsText, descriptionText);
            string customMonosaccharidePath = GlobalVariables.CustomMonosaccharidePath;
            try
            {
                Directory.CreateDirectory(Path.GetDirectoryName(customMonosaccharidePath));
                if (!File.Exists(customMonosaccharidePath))
                {
                    File.WriteAllLines(customMonosaccharidePath, new[] { MonoSaccharidesHeader, line });
                }
                else
                {
                    // AppendAllLines never checks whether the file already ends in a newline --
                    // if a hand-edited TSV's last byte isn't one (easy to end up with outside the
                    // shipped file, which is exactly who this window is for), the new row gets
                    // glued onto the end of the last existing line instead of starting its own.
                    string existing = File.ReadAllText(customMonosaccharidePath);
                    if(existing.Length > 0 && !existing.EndsWith("\n"))
                    {
                        File.AppendAllText(customMonosaccharidePath, Environment.NewLine);
                    }
                    File.AppendAllLines(customMonosaccharidePath, new[] { line });
                }
                return null;
            }
            catch (Exception ex)
            {
                return $"The monosaccharide is available for this session, but could not be saved to file for future sessions: {ex.Message}";
            }
        }


        /// <summary>
        /// The two ways a glycan database line can be written. A database file is read entirely as one or
        /// the other -- the format is inferred from its first data line and never declared -- so an entry
        /// being added has to agree with whatever is already in the file.
        /// </summary>
        public enum GlycanLineFormat
        {
            /// <summary>Nested single-character codes, e.g. <c>(N(H(A)))</c>. What OGlycan.gdb uses.</summary>
            Structure,

            /// <summary>Name-and-count, e.g. <c>HexNAc(2)Hex(5)</c>. What NGlycan.gdb uses.</summary>
            Composition
        }

        /// <summary>
        /// Validates one glycan the user has typed and appends it to their custom glycan database, in the
        /// format that database is already written in.
        ///
        /// <para>
        /// Return contract, and it differs deliberately from <see cref="PersistCustomMonosaccharide"/>:
        /// this method throws <see cref="MetaMorpheusException"/> for every failure and returns normally
        /// only once the glycan is on disk. There is no registered-but-not-saved middle state to report,
        /// because a glycan database is not held in memory the way a monosaccharide is -- the search reads
        /// the file when the engine is built. Which is also why an entry added now is picked up by a search
        /// run in this same session; only MetaDraw's glycan list waits for a restart.
        /// </para>
        /// </summary>
        /// <param name="glycanText">The glycan, as typed: a structure or a composition.</param>
        /// <param name="databasePath">The custom database to append to. Created if it is not there.</param>
        /// <param name="isOGlycan">
        /// Whether this is an O-glycan database, so the entry is validated the same way the search will
        /// read it and a glycan that parses here cannot fail to parse there.
        /// </param>
        public static void PersistCustomGlycan(string glycanText, string databasePath, bool isOGlycan)
        {
            string entry = (glycanText ?? string.Empty).Trim();
            string fileName = Path.GetFileName(databasePath);

            if (entry.Length == 0)
            {
                throw new MetaMorpheusException("Could not add the glycan: no glycan was given.");
            }
            if (IsCommentOrBlank(entry))
            {
                throw new MetaMorpheusException(
                    $"Could not add the glycan: \"{entry}\" is a comment, not a glycan. Lines beginning with '#' are ignored when the database is read.");
            }
            if (entry.IndexOfAny(new[] { '\t', '\n', '\r' }) >= 0)
            {
                throw new MetaMorpheusException(
                    "Could not add the glycan: one glycan per line, with no tabs or line breaks in it.");
            }

            GlycanLineFormat format = ValidateGlycanLine(entry, fileName, isOGlycan);

            // The format is a property of the whole file, so an entry that disagreed with what is already
            // there would silently change how every OTHER line in the file is read.
            GlycanLineFormat? existingFormat = FormatOfExistingEntries(databasePath);
            if (existingFormat.HasValue && existingFormat.Value != format)
            {
                throw new MetaMorpheusException(
                    $"Could not add the glycan to '{fileName}': that file is written in {existingFormat.Value.ToString().ToLowerInvariant()} format " +
                    $"and \"{entry}\" is {format.ToString().ToLowerInvariant()} format. A glycan database is read entirely as one format or the other, " +
                    "so the two cannot be mixed in one file. Convert the entry, or keep it in a database of its own.");
            }

            string duplicate = ExistingEntryFor(databasePath, entry, format);
            if (duplicate != null)
            {
                string sameGlycan = duplicate == entry ? "" : $", which is the same glycan as \"{entry}\"";
                throw new MetaMorpheusException($"Could not add the glycan: '{fileName}' already contains \"{duplicate}\"{sameGlycan}.");
            }

            try
            {
                string directory = Path.GetDirectoryName(databasePath);
                if (!string.IsNullOrEmpty(directory))
                {
                    Directory.CreateDirectory(directory);
                }

                if (!File.Exists(databasePath))
                {
                    File.WriteAllLines(databasePath, new[] { entry });
                }
                else
                {
                    // AppendAllLines never checks whether the file already ends in a newline. A hand-edited
                    // .gdb very often does not -- NGlycan_ForNoSearch.gdb in the test data does not -- and
                    // without this the new glycan is glued onto the end of the last one, corrupting both.
                    string existing = File.ReadAllText(databasePath);
                    if (existing.Length > 0 && !existing.EndsWith("\n"))
                    {
                        File.AppendAllText(databasePath, Environment.NewLine);
                    }
                    File.AppendAllLines(databasePath, new[] { entry });
                }
            }
            catch (Exception ex)
            {
                throw new MetaMorpheusException($"Could not save the glycan to '{databasePath}': {ex.Message}", ex);
            }
        }

        /// <summary>
        /// Works out which format a glycan line is written in, and checks that the parser the search will
        /// use actually accepts it -- so a glycan the user is told was added cannot fail to load later.
        /// </summary>
        /// <exception cref="MetaMorpheusException">The line is neither format, or does not parse.</exception>
        public static GlycanLineFormat ValidateGlycanLine(string entry, string fileName, bool isOGlycan)
        {
            if (entry.StartsWith("(", StringComparison.Ordinal))
            {
                ValidateStructureEntry(entry, isOGlycan);
                return GlycanLineFormat.Structure;
            }

            ValidateCompositionEntry(entry);
            return GlycanLineFormat.Composition;
        }

        /// <summary>
        /// The most monosaccharides a glycan entered through <see cref="PersistCustomGlycan"/> may hold --
        /// twice the largest glycan in any database MetaMorpheus ships (20, in NGlycan.gdb).
        /// </summary>
        /// <remarks>
        /// Zero is refused because it means nothing; this is the other end of the same range. The child ions
        /// of a glycan grow combinatorially with its size, and they are built when a search first reads the
        /// database: HexNAc(255)Hex(255)NeuAc(255) takes nearly four minutes there, with no sign of why. The
        /// loader applies no cap -- a database that already loads is not refused after the fact.
        /// </remarks>
        public const int MaxMonosaccharidesPerGlycan = 40;

        private static void ValidateStructureEntry(string entry, bool isOGlycan)
        {
            // The same checks the loader runs -- characters, balance, one root -- so a structure is refused
            // here, naming what was typed, rather than being told it was added and failing at search time.
            string problem = StructureProblem(entry);
            if (problem != null)
            {
                throw new MetaMorpheusException($"Could not add the glycan \"{entry}\": {problem}");
            }

            // Counted before Struct2Glycan, which is where the combinatorial cost is paid.
            int residues = entry.Count(c => c != '(' && c != ')');
            if (residues > MaxMonosaccharidesPerGlycan)
            {
                throw new MetaMorpheusException(TooLarge(entry, residues));
            }

            try
            {
                // Parsed exactly as LoadStructureGlycan would parse it: Struct2Glycan is the only thing that
                // knows whether the nesting describes a tree it can actually build.
                List<Glycan> parsed = Glycan.Struct2Glycan(entry, 1, isOGlycan);

                // "()" parses happily into a glycan of nothing, whose mass is zero. Searching for it is
                // meaningless and it would widen every box it landed in, so it is refused here rather than
                // left for the user to wonder about. The empty case is folded in: Struct2Glycan either
                // throws or yields glycans, so "nothing came back" and "nothing in what came back" are the
                // same answer to the user.
                if (parsed == null || parsed.Count == 0 || parsed.All(g => g.Kind.Sum(count => (int)count) == 0))
                {
                    throw new MetaMorpheusException(
                        $"Could not add the glycan \"{entry}\": it contains no monosaccharides.");
                }
            }
            catch (MetaMorpheusException)
            {
                throw;
            }
            catch (Exception ex)
            {
                throw new MetaMorpheusException($"Could not add the glycan \"{entry}\": {ex.Message}", ex);
            }
        }

        /// <summary>
        /// Why a structure line cannot be read, or null when it can. Only the reason: the loader and the
        /// entry validator each say where it came from -- a file and a line, or what the user typed.
        /// </summary>
        private static string StructureProblem(string structure)
        {
            // Struct2Glycan silently miscounts a character it has no mass for, so an unknown monosaccharide
            // code would be searched as a lighter glycan than the one written. Consults the live registry,
            // so codes declared in MonosaccharidesCustom.tsv are accepted too.
            foreach (char c in structure)
            {
                if (c != '(' && c != ')' && !Glycan.CharMassDic.ContainsKey(c))
                {
                    return $"Unrecognized character '{c}'. Allowed: parentheses and one of {string.Concat(Glycan.CharMassDic.Keys)}. " +
                        "A monosaccharide MetaMorpheus does not ship with must be declared in MonosaccharidesCustom.tsv first.";
                }
            }

            int depth = 0;
            int firstTreeEnds = -1;
            for (int i = 0; i < structure.Length; i++)
            {
                if (structure[i] == '(')
                {
                    depth++;
                }
                else if (structure[i] == ')')
                {
                    depth--;
                    if (depth < 0)
                    {
                        return "The parentheses do not balance -- a ')' closes a branch that was never opened.";
                    }
                    if (depth == 0 && firstTreeEnds < 0)
                    {
                        firstTreeEnds = i;
                    }
                }
            }
            if (depth != 0)
            {
                return $"The parentheses do not balance -- {depth} branch(es) are left open.";
            }

            // Struct2Node reads one tree and stops, so "(N)(H)" -- which balances -- would be searched as
            // (N) with the H quietly gone. A second root is a second tree.
            if (firstTreeEnds >= 0 && firstTreeEnds < structure.Length - 1)
            {
                return "A structure is a single tree with one root, and this has more than one: " +
                    $"\"{structure.Substring(0, firstTreeEnds + 1)}\" is complete before \"{structure.Substring(firstTreeEnds + 1)}\" begins.";
            }

            return null;
        }

        /// <summary>
        /// A composition: a monosaccharide name and a count, repeated, e.g. HexNAc(2)Hex(5). Shared so the
        /// loader and the entry validator cannot drift apart about what one looks like -- which is exactly
        /// how the format-detection bug got in. The name is anything but a parenthesis, because that is all
        /// MonosaccharidesCustom.tsv asks of one: Hex-6P and Neu5,9Ac2 are legal names.
        /// </summary>
        private const string CompositionShape = @"^(?:(?<name>[^()]+)\((?<count>\d+)\))+$";

        /// <summary>
        /// The part of a line in a composition database that is the glycan: before the first tab (the
        /// shipped .txt databases carry name and mass columns there), before any '#' note, and up to the
        /// last ')' -- so a column lined up with spaces is dropped the way a tabbed one is. String2Kind
        /// always ignored whatever followed the last ')'; this keeps that tolerance explicit.
        /// </summary>
        private static string CompositionPart(string line)
        {
            string glycan = line.Split('\t')[0];
            int note = glycan.IndexOf('#');
            if (note >= 0)
            {
                glycan = glycan.Substring(0, note);
            }
            int lastClose = glycan.LastIndexOf(')');
            if (lastClose >= 0)
            {
                glycan = glycan.Substring(0, lastClose + 1);
            }
            return glycan.Trim();
        }

        /// <summary>
        /// Parse a composition into its kind[], or say why it cannot be. Throws FormatException carrying only
        /// the reason, for the caller to say where the composition came from.
        /// </summary>
        /// <remarks>
        /// Unlike String2Kind, which throws KeyNotFoundException on an unknown name, quietly accepts a
        /// trailing unclosed group, and lets a second name for the same monosaccharide (Fuc and dHex are one
        /// slot) overwrite the first.
        /// </remarks>
        private static byte[] ParseComposition(string composition)
        {
            Match match = Regex.Match(composition, CompositionShape);
            if (!match.Success)
            {
                throw new FormatException(
                    "It is neither a structure -- which starts with '(', e.g. (N(H(A))) -- " +
                    "nor a composition, which is a name and a count repeated, e.g. HexNAc(2)Hex(5).");
            }

            CaptureCollection names = match.Groups["name"].Captures;
            CaptureCollection counts = match.Groups["count"].Captures;
            byte[] kind = new byte[Glycan.KindCapacity];
            var nameInSlot = new Dictionary<int, string>();

            for (int i = 0; i < names.Count; i++)
            {
                string name = names[i].Value;
                if (!Glycan.NameCharDic.TryGetValue(name, out var code))
                {
                    throw new FormatException(
                        $"'{name}' is not a monosaccharide MetaMorpheus knows. Known: {string.Join(", ", Glycan.NameCharDic.Keys)}. " +
                        "A monosaccharide it does not ship with must be declared in MonosaccharidesCustom.tsv first.");
                }
                if (nameInSlot.TryGetValue(code.Item2, out string earlier))
                {
                    throw new FormatException(earlier == name
                        ? $"'{name}' appears more than once. Give each monosaccharide a single total count."
                        : $"'{earlier}' and '{name}' are two names for the same monosaccharide. Give it once, with a single total count.");
                }
                if (!byte.TryParse(counts[i].Value, NumberStyles.Integer, CultureInfo.InvariantCulture, out byte count))
                {
                    throw new FormatException($"The count for '{name}' must be a whole number between 0 and 255.");
                }
                nameInSlot[code.Item2] = name;
                kind[code.Item2] = count;
            }

            return kind;
        }

        private static void ValidateCompositionEntry(string entry)
        {
            byte[] kind;
            try
            {
                kind = ParseComposition(entry);
            }
            catch (FormatException ex)
            {
                throw new MetaMorpheusException($"Could not add the glycan \"{entry}\": {ex.Message}");
            }

            int total = kind.Sum(count => (int)count);

            // Same reason as the structure case: a composition that totals nothing is a zero-mass glycan.
            if (total == 0)
            {
                throw new MetaMorpheusException(
                    $"Could not add the glycan \"{entry}\": it contains no monosaccharides.");
            }
            if (total > MaxMonosaccharidesPerGlycan)
            {
                throw new MetaMorpheusException(TooLarge(entry, total));
            }
        }

        private static string TooLarge(string entry, int residues) =>
            $"Could not add the glycan \"{entry}\": it has {residues} monosaccharides, and a glycan may have at most " +
            $"{MaxMonosaccharidesPerGlycan}. The ions a glycan search builds grow combinatorially with its size, so one this " +
            "large would stall the first search that reads the database.";

        /// <summary>
        /// The format the entries already in a database are written in, or null when it holds no entries
        /// yet -- a file that is missing, or one that is all banner, which is what a freshly seeded custom
        /// database is.
        /// </summary>
        private static GlycanLineFormat? FormatOfExistingEntries(string databasePath)
        {
            if (!File.Exists(databasePath))
            {
                return null;
            }

            foreach (string line in File.ReadLines(databasePath))
            {
                if (IsCommentOrBlank(line))
                {
                    continue;
                }
                return line.TrimStart().StartsWith("(", StringComparison.Ordinal)
                    ? GlycanLineFormat.Structure
                    : GlycanLineFormat.Composition;
            }

            return null;
        }

        /// <summary>
        /// The line already in the database that is the same glycan as the entry, or null. A composition is a
        /// set of counts, so HexNAc(1)Hex(1), Hex(1)HexNAc(1), HexNAc(1)Hex(1)Fuc(0) and -- Fuc and dHex being
        /// one slot -- Fuc(1) and dHex(1) are compared by what they parse to. A structure is compared as
        /// written: two trees of one composition are different glycans.
        /// </summary>
        private static string ExistingEntryFor(string databasePath, string entry, GlycanLineFormat format)
        {
            if (!File.Exists(databasePath))
            {
                return null;
            }

            byte[] entryKind = format == GlycanLineFormat.Composition ? ParseComposition(entry) : null;

            foreach (string line in File.ReadLines(databasePath))
            {
                if (IsCommentOrBlank(line))
                {
                    continue;
                }

                if (entryKind == null)
                {
                    // Split on tab first: a database may carry columns after the glycan.
                    string structure = line.Split('\t')[0].Trim();
                    if (structure.Equals(entry, StringComparison.Ordinal))
                    {
                        return structure;
                    }
                    continue;
                }

                string composition = CompositionPart(line);
                try
                {
                    if (ParseComposition(composition).SequenceEqual(entryKind))
                    {
                        return composition;
                    }
                }
                catch (FormatException)
                {
                    // A line this cannot read is the loader's to report, not a reason to refuse the entry.
                }
            }

            return null;
        }

        /// <summary>
        /// Ensure the MonosaccharidesCustom.tsv exists in the directory. If the file is missing, 
        /// write the embedded fully documented template—instructions, column spec, the built-in name/code table,
        /// worked examples — with the header row as its single non-comment line; do nothing if it already exists.
        /// </summary>
        /// <param name="path">
        /// The destination path — normally GlobalVariables.CustomMonosaccharidePath.
        /// </param>
        public static void EnsureCustomMonosaccharideFileExists(string path) 
        {
            if (!File.Exists(path)) 
            { 
                try
                {
                    // Make sure the directory exists before writing the file, however, DataDir is created by
                    // SetUpDataDirectory before this runs; this is defensive only.
                    Directory.CreateDirectory(Path.GetDirectoryName(path));

                    // Non-installer/portable runs (or a user-specified DataDir) may still have the old
                    // Glycan_Mods\MonosaccharidesCustom.tsv from before this fix. Carry it over once so
                    // those users don't silently lose custom entries; leave the old file alone.
                    string legacyPath = Path.Combine(Path.GetDirectoryName(path), "Glycan_Mods", "MonosaccharidesCustom.tsv");
                    if (File.Exists(legacyPath))
                    {
                        File.Copy(legacyPath, path);
                        return;
                    }

                    // The default template—instructions is embedded in the DLL so it survives
                    // install/repair/upgrade regardless of what the installer does to Glycan_Mods --
                    // same pattern as RnaMods.txt (GlobalVariables.LoadRnaModifications) and the
                    // mzLib-embedded default protease/rnase templates (GlobalVariables.LoadDigestionAgents).
                    var assembly = typeof(GlycanDatabase).Assembly;
                    using var stream = assembly.GetManifestResourceStream("EngineLayer.Glycan_Mods.MonosaccharidesCustom.tsv");
                    using var reader = new StreamReader(stream);
                    File.WriteAllText(path, reader.ReadToEnd());
                }
                catch (Exception ex)
                {
                    throw new MetaMorpheusException($"Could not create the custom monosaccharide file '{path}': {ex.Message}", ex);
                }
            }
        }


        /// <summary>
        /// Load composition format Glycan database, then convert to kind format followed by generating the glycan object.
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="ToGenerateIons"></param>
        /// <param name="IsOGlycanSearch"></param>
        /// <returns>The glycan collection </returns>
        public static IEnumerable<Glycan> LoadKindGlycan(string filePath, bool ToGenerateIons, bool IsOGlycan)
        {
            using (StreamReader lines = new StreamReader(filePath))
            {
                int id = 1;
                int lineNumber = 0;
                while (lines.Peek() != -1)
                {
                    string rawLine = lines.ReadLine();
                    lineNumber++;

                    // Skipped explicitly rather than left to the shape test below: a comment that quotes a
                    // composition -- which the documented template does, to show the format -- would otherwise
                    // reach String2Kind and die on a dictionary lookup naming neither the file nor the line.
                    if (IsCommentOrBlank(rawLine))
                    {
                        continue;
                    }

                    string line = CompositionPart(rawLine);

                    // A line with no parenthesis in it has not begun to be a composition -- a column header,
                    // or stray text -- and is skipped, as it always was. The test used to be "does it contain
                    // Hex or HexNAc", which also threw away, without a word, every composition built from
                    // neither: NeuAc(2)Fuc(1) is one, and the entry validator accepts it.
                    if (line.IndexOf('(') < 0)
                    {
                        continue;
                    }

                    byte[] kind;
                    try
                    {
                        kind = ParseComposition(line);  // Convert the database string to kind[] format (byte array).
                    }
                    catch (FormatException ex)
                    {
                        // A line that has started to be a composition and is not one is named rather than
                        // dropped: skipping it would search a database the user did not write, and say
                        // nothing. Named the way a bad structure line is.
                        throw new MetaMorpheusException(
                            $"Could not parse glycan composition in '{Path.GetFileName(filePath)}' at line {lineNumber}: \"{line}\". {ex.Message}",
                            ex);
                    }

                    if (IsOGlycan) // Load the oGlycan with two different motifs : S and T
                    {
                        var oGlycan_S = new Glycan(kind, "S", GlycanType.O_glycan); // Use the kind[] to create a glycan object.  
                        oGlycan_S.GlyId = id;
                        id++;
                        if (ToGenerateIons)
                        {
                            oGlycan_S.Ions = OGlycanCompositionCombinationChildIons(kind);
                        }
                        yield return oGlycan_S; // Output the first glycan  

                        var oGlycan_T = new Glycan(kind, "T", GlycanType.O_glycan); // Use the kind[] to create a glycan object.  
                        oGlycan_T.GlyId = id;
                        id++;
                        if (ToGenerateIons)
                        {
                            oGlycan_T.Ions = OGlycanCompositionCombinationChildIons(kind);
                        }
                        yield return oGlycan_T; // Output the second glycan  
                    }
                    else // Load the N-glycan with one motif : N
                    {
                        var nGlycan_Nxs = new Glycan(kind, "Nxs", GlycanType.N_glycan); // Use the kind[] to create a glycan object.
                        nGlycan_Nxs.GlyId = id;
                        id++;
                        if (ToGenerateIons)
                        {
                            nGlycan_Nxs.Ions = OGlycanCompositionCombinationChildIons(kind);
                        }
                        yield return nGlycan_Nxs;


                        var nGlycan_Nxt = new Glycan(kind, "Nxt", GlycanType.N_glycan); // Use the kind[] to create a glycan object.
                        nGlycan_Nxt.GlyId = id;
                        id++;
                        if (ToGenerateIons)
                        {
                            nGlycan_Nxt.Ions = OGlycanCompositionCombinationChildIons(kind);
                        }
                        yield return nGlycan_Nxt;
                    }
                }
            }
        }

        /// <summary>
        /// Convert the glycan string to Kind array
        /// </summary>
        /// <param name="line"> ex. HexNAc(2)Hex(5)NeuAc(1)Fuc(1) </param>
        /// <returns> The glycan Kind List ex. [2, 5, 0, 0, 1, 0, 0, 0, 0, 1] </returns>
        public static byte[] String2Kind(string line) 
        {
            byte[] kind = new byte[Glycan.KindCapacity];
            var x = line.Split(new char[] { '(', ')' });
            int i = 0;
            while (i < x.Length - 1)
            {
                kind[Glycan.NameCharDic[x[i]].Item2] = byte.Parse(x[i + 1]);
                i = i + 2;
            }

            return kind;
        }

        /// <summary>
        /// Load structured format Glycan database and generate the glycan object.
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="IsOGlycan"></param>
        /// <returns> The Glycan object collection </returns>
        public static IEnumerable<Glycan> LoadStructureGlycan(string filePath, bool IsOGlycan)
        {
            using (StreamReader glycans = new StreamReader(filePath))
            {
                int id = 1;
                int lineNumber = 0;
                while (glycans.Peek() != -1)
                {
                    string line = glycans.ReadLine();   // Read the line from the database file. Ex. (N(H(A))(A))
                    lineNumber++;

                    // A '#' banner or a blank spacer is not a glycan. Without this, the seeded template --
                    // and any database a user has annotated -- reaches the structure check and throws on
                    // the '#' itself during startup, before any window opens.
                    if (IsCommentOrBlank(line))
                    {
                        continue;
                    }

                    // Characters the parser would silently miscount as zero mass, parentheses that do not
                    // balance, and a second root Struct2Node would silently drop.
                    string structureProblem = StructureProblem(line.Trim());
                    if (structureProblem != null)
                    {
                        throw new MetaMorpheusException(
                            $"Could not parse glycan structure in '{Path.GetFileName(filePath)}' at line {lineNumber}: \"{line.Trim()}\". {structureProblem}");
                    }

                    // For each glycan, two versions will be generated:
                    // For O-glycan, one modified on serine (S), and the other on threonine (T).
                    // For N-glycan, one modified on N-glycosylation on motif Asn-X-Ser(Nxs), and the other on Asn-X-Thr(Nxt).
                    foreach (var glycan in Glycan.Struct2Glycan(line, id, IsOGlycan)) // Modify the line to handle multiple Glycan objects returned by Struct2Glycan.
                    {
                        yield return glycan;
                    }
                    id = id + 2; // Each line will generate two glycan objects
                }
            }
        }

        // A line is treated as a comment if its first non-whitespace character is '#'.
        // Blank/whitespace-only lines are also skipped so users can space out their database files.
        private static bool IsCommentOrBlank(string line)
        {
            if (string.IsNullOrWhiteSpace(line))
                return true;
            // Microsoft Excel sometimes wraps lines in double-quotes when saving as TSV/CSV.
            // A comment line that originally starts with '#' may appear as "# MY COMMENT after
            // a round-trip through Excel.  Strip a leading '"' before checking for '#'.
            string trimmed = line.TrimStart().TrimStart('"');
            return trimmed.StartsWith("#");
        }

        //This function build fragments based on the general core of NGlyco fragments. 
        //From https://github.com/mobiusklein/glycopeptidepy/structure/fragmentation_strategy/glycan.py#L408
        //The fragment generation is not as good as structure based method. So it is better to use a structure based N-Glycan database.
        // The function is used to load the database from the different formats, but we don't use it now.
        public static List<GlycanIon> NGlycanCompositionFragments(byte[] kind, bool isfucExtended = false)
        {
            int glycan_mass = Glycan.GetMass(kind);

            // int core_count = 1;
            int iteration_count = 0;
            int hexnac_Core = 2;
            int hexose_Core = 3;
            bool extended = true;
            bool extended_fucosylation = isfucExtended;

            int fuc_count = kind[4];
            int xyl_count = kind[9];
            int hexnac_total = kind[1];
            int hexose_total = kind[0];

            List<GlycanIon> glycanIons = new List<GlycanIon>();

            int base_hexnac = Math.Min(hexnac_total, hexnac_Core); // base_hexnac is the first priority hexnac count, they all come from the core.
            for (int hexnac_count = 0; hexnac_count < base_hexnac + 1 ; hexnac_count++)
            {
                if (hexnac_count == 0)
                {
                    byte[] startKind = new byte[Glycan.KindCapacity];
                    startKind[1] = (byte)hexnac_count;
                    string glycanName = Glycan.GetKindString(startKind);
                    GlycanIon glycanIon = new GlycanIon(glycanName, 8303819, startKind, glycan_mass - 8303819);
                    glycanIons.Add(glycanIon);
                }
                else if (hexnac_count == 1)
                {
                    GlycanIon glycanIon = GenerateGlycanIon(0, (byte)hexnac_count, 0, 0, glycan_mass);

                    glycanIons.Add(glycanIon);

                    if (iteration_count < fuc_count)
                    {
                        GlycanIon fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);

                        glycanIons.Add(fuc_glycanIon);
                    }
                }
                else if (hexnac_count == 2)
                {
                    GlycanIon glycanIon = GenerateGlycanIon(0, (byte)hexnac_count, 0, 0, glycan_mass);
                    glycanIons.Add(glycanIon);

                    if (!extended_fucosylation)
                    {
                        if (iteration_count < fuc_count)
                        {
                            GlycanIon fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);
                            glycanIons.Add(fuc_glycanIon);

                            if (iteration_count < xyl_count)
                            {
                                GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                glycanIons.Add(xyl_fuc_glycanIon);
                            }
                        }
                    }
                    else if (fuc_count > 0)
                    {
                        GlycanIon fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);
                        glycanIons.Add(fuc_glycanIon);

                        for (int add_fuc_count = 2; add_fuc_count <= fuc_count; add_fuc_count++)
                        {
                            GlycanIon add_fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);
                            glycanIons.Add(add_fuc_glycanIon);
                        }

                        if (iteration_count < xyl_count)
                        {
                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                            glycanIons.Add(xyl_fuc_glycanIon);
                        }
                    }

                    if (iteration_count < xyl_count)
                    {
                        GlycanIon xyl_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 0, 1, glycan_mass);
                        glycanIons.Add(xyl_glycanIon);
                    }


                    int base_hexose = Math.Min(hexose_total, hexose_Core); // base_hexose is the first priority hexose count, they all come from the core.
                    for (int hexose_count = 1; hexose_count <= base_hexose + 1; hexose_count++)
                    {
                        GlycanIon hexose_glycanIon = GenerateGlycanIon((byte)hexose_count, (byte)hexnac_count, 0, 0, glycan_mass);
                        glycanIons.Add(hexose_glycanIon);

                        if (!extended_fucosylation)
                        {
                            if (iteration_count < fuc_count)
                            {
                                GlycanIon fuc_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 1, 0, glycan_mass);
                                glycanIons.Add(fuc_glycanIon);

                                if (iteration_count < xyl_count)
                                {
                                    GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                    glycanIons.Add(xyl_fuc_glycanIon);
                                }
                            }                           
                        }
                        else if (fuc_count > 0)
                        {
                            GlycanIon fuc_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 1, 0, glycan_mass);
                            glycanIons.Add(fuc_glycanIon);

                            for (int add_fuc_count = 2; add_fuc_count <= fuc_count; add_fuc_count++)
                            {
                                GlycanIon add_fuc_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 1, 0, glycan_mass);
                                glycanIons.Add(add_fuc_glycanIon);
                            }

                            if (iteration_count < xyl_count)
                            {
                                GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                glycanIons.Add(xyl_fuc_glycanIon);
                            }
                        }

                        if (iteration_count < xyl_count)
                        {
                            GlycanIon xyl_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 0, 1, glycan_mass);
                            glycanIons.Add(xyl_glycanIon);
                        }

                        if (hexose_count == hexose_Core && hexnac_count >= hexnac_Core  && extended) //After the core motif has been exhausted, speculatively add on the remaining core monosaccharides sequentially until exhausted.
                        {
                            for (int extra_hexnac_count = 0; extra_hexnac_count < hexnac_total - hexnac_count + 1; extra_hexnac_count++)
                            {
                                if (extra_hexnac_count + hexnac_count > hexnac_total) // this part is doesn't make sense, because the hexnac_count cannot be larger than total-hexnac
                                {
                                    continue;
                                }

                                if (extra_hexnac_count > 0)
                                {
                                    GlycanIon new_glycanIon = GenerateGlycanIon((byte)hexose_count, (byte)(hexnac_count + extra_hexnac_count), 0, 0, glycan_mass);

                                    glycanIons.Add(new_glycanIon);

                                    if (!extended_fucosylation)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                        glycanIons.Add(fuc_glycanIon);

                                        if (iteration_count < xyl_count)
                                        {
                                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                            glycanIons.Add(xyl_fuc_glycanIon);
                                        }
                                    }
                                    else if (fuc_count > 0)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                        glycanIons.Add(fuc_glycanIon);

                                        for (int add_fuc_count = 2; add_fuc_count <= fuc_count; add_fuc_count++)
                                        {
                                            GlycanIon add_fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                            glycanIons.Add(add_fuc_glycanIon);
                                        }

                                        if (iteration_count < xyl_count)
                                        {
                                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                            glycanIons.Add(xyl_fuc_glycanIon);
                                        }
                                    }

                                    if (iteration_count < xyl_count)
                                    {
                                        GlycanIon xyl_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 0, 1, glycan_mass);
                                        glycanIons.Add(xyl_glycanIon);
                                    }

                                }

                                for (int extra_hexose_count = 1; extra_hexose_count < hexose_total - hexose_Core + 1; extra_hexose_count++)
                                {
                                    if (extra_hexose_count + hexose_count > hexose_total) // this part is doesn't make sense, because the hexnac_count cannot be larger than total-hexnac
                                    {
                                        continue;
                                    }

                                    GlycanIon new_glycanIon = GenerateGlycanIon((byte)(hexose_count + extra_hexose_count), (byte)(hexnac_count + extra_hexnac_count), 0, 0, glycan_mass);

                                    glycanIons.Add(new_glycanIon);

                                    if (!extended_fucosylation)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                        glycanIons.Add(fuc_glycanIon);

                                        if (iteration_count < xyl_count)
                                        {
                                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                            glycanIons.Add(xyl_fuc_glycanIon);
                                        }
                                    }
                                    else if (fuc_count > 0)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                        glycanIons.Add(fuc_glycanIon);

                                        for (int add_fuc_count = 2; add_fuc_count <= fuc_count; add_fuc_count++)
                                        {
                                            GlycanIon add_fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);
                                            glycanIons.Add(add_fuc_glycanIon);
                                        }

                                        if (iteration_count < xyl_count)
                                        {
                                            GlycanIon xyl_fuc_glycanIon = ExtendGlycanIon(fuc_glycanIon, 0, 0, 0, 1, glycan_mass);
                                            glycanIons.Add(xyl_fuc_glycanIon);
                                        }
                                    }

                                    if (iteration_count < xyl_count)
                                    {
                                        GlycanIon xyl_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 0, 1, glycan_mass);
                                        glycanIons.Add(xyl_glycanIon);
                                    }

                                }
                            }
                        }
                    }

                }


            }

            return glycanIons;
        }

        private static GlycanIon GenerateGlycanIon(byte hexose_count, byte hexnac_count, byte fuc_count, byte xyl_count, int glycan_mass)
        {
            byte[] ionKind = new byte[Glycan.KindCapacity];
            ionKind[0] = hexose_count;
            ionKind[1] = hexnac_count;
            ionKind[4] = fuc_count;
            ionKind[9] = xyl_count;

            int ionMass = Glycan.GetMass(ionKind);

            String glycanName = Glycan.GetKindString(ionKind);

            GlycanIon glycanIon = new GlycanIon(glycanName, ionMass, ionKind, glycan_mass - ionMass);

            return glycanIon;
        }

        private static GlycanIon ExtendGlycanIon(GlycanIon glycanIon, byte hexose_count, byte hexnac_count, byte fuc_count, byte xyl_count, int glycan_mass)
        {
            byte[] ionKind = glycanIon.IonKind;
            ionKind[0] += hexose_count;
            ionKind[1] += hexnac_count;
            ionKind[4] += fuc_count;
            ionKind[9] += xyl_count;

            int ionMass = Glycan.GetMass(ionKind);
            string glycanName = Glycan.GetKindString(ionKind);

            GlycanIon extend_glycanIon = new GlycanIon(glycanName, ionMass, ionKind, glycan_mass - ionMass);

            return extend_glycanIon;
        }

        //This function build fragments based on the general core of OGlyco fragments. 
        //From https://github.com/mobiusklein/glycopeptidepy/structure/fragmentation_strategy/glycan.py
        //The fragment generation is not as good as structure based method. So it is better to use a structure based O-Glycan database.
        // We don't use this function now, alternatively, we use the 'OGlycanCompositionCombinationChildIons'.
        public static List<GlycanIon> OGlycanCompositionFragments(byte[] kind)
        {
            List<GlycanIon> glycanIons = new List<GlycanIon>();

            int glycan_mass = Glycan.GetMass(kind);

            int iteration_count = 0;
            bool extended = true;

            int fuc_count = kind[4];
            int hexnac_total = kind[1];
            int hexose_total = kind[0];

            for (int hexnac_count = 0; hexnac_count < 3; hexnac_count++)
            {
                if (hexnac_total < hexnac_count)
                {
                    continue;
                }


                if (hexnac_count >= 1)
                {
                    GlycanIon glycanIon = GenerateGlycanIon(0, (byte)hexnac_count, 0, 0, glycan_mass);

                    glycanIons.Add(glycanIon);

                    if (iteration_count < fuc_count)
                    {
                        GlycanIon fuc_glycanIon = ExtendGlycanIon(glycanIon, 0, 0, 1, 0, glycan_mass);

                        glycanIons.Add(fuc_glycanIon);
                    }

                    for (int hexose_count = 0; hexose_count < 2; hexose_count++)
                    {
                        if (hexose_total < hexose_count)
                        {
                            continue;
                        }

                        if (hexose_count > 0)
                        {
                            GlycanIon hexose_glycanIon = GenerateGlycanIon((byte)hexose_count, (byte)hexnac_count, 0, 0, glycan_mass);
                            glycanIons.Add(hexose_glycanIon);

                            if (iteration_count < fuc_count)
                            {
                                GlycanIon fuc_glycanIon = ExtendGlycanIon(hexose_glycanIon, 0, 0, 1, 0, glycan_mass);

                                glycanIons.Add(fuc_glycanIon);
                            }
                        }

                        // After the core motif has been exhausted, speculatively add on the remaining core monosaccharides sequentially until exhausted.

                        if (extended && hexnac_total - hexnac_count >= 0)
                        {
                            for (int extra_hexnac_count = 0; extra_hexnac_count  < hexnac_total - hexnac_count + 1; extra_hexnac_count ++)
                            {
                                if (extra_hexnac_count > 0)
                                {
                                    GlycanIon new_glycanIon = GenerateGlycanIon((byte)hexose_count, (byte)(hexnac_count + extra_hexnac_count), 0, 0, glycan_mass);

                                    glycanIons.Add(new_glycanIon);


                                    if (iteration_count < fuc_count)
                                    {
                                        GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);

                                        glycanIons.Add(fuc_glycanIon);
                                    }

                                }

                                if (hexose_total > hexose_count && hexose_count > 0)
                                {
                                    for (int extra_hexose_count = 0; extra_hexose_count < hexose_total - hexose_count; extra_hexose_count++)
                                    {
                                        if (extra_hexose_count > 0 && extra_hexose_count + hexose_count >0)
                                        {

                                            GlycanIon new_glycanIon = GenerateGlycanIon((byte)(hexose_count + extra_hexose_count), (byte)(hexnac_count + extra_hexnac_count), 0, 0, glycan_mass);

                                            glycanIons.Add(new_glycanIon);


                                            if (iteration_count < fuc_count)
                                            {
                                                GlycanIon fuc_glycanIon = ExtendGlycanIon(new_glycanIon, 0, 0, 1, 0, glycan_mass);

                                                glycanIons.Add(fuc_glycanIon);
                                            }
                                        }
                                    }
                                }
                            }
                        }

                    }
                }

            }


            return glycanIons;
        }

        /// <summary>
        /// Generate some child ions based on the kind array. The kind array is the combination of the monosaccharides then filter by the rules.
        /// </summary>
        /// <param name="kind"> glycan Kind[]</param>
        /// <returns> The glycanIon collection </returns>
        public static List<GlycanIon> OGlycanCompositionCombinationChildIons(byte[] kind)
        {
            List<GlycanIon> glycanIons = new List<GlycanIon>();

            int glycan_mass = Glycan.GetMass(kind);

            List<byte[]> _kinds = new List<byte[]>();
            HashSet<string> _keys = new HashSet<string>();

            _kinds.Add((byte[])kind.Clone());
            _GetCombinations(kind, _kinds, _keys);

            foreach (var k in _kinds)
            {
                //Rules to build OGlycan child ions. Filter the kind array which doesn't meet the rules.
                //At least one HexNAc
                if (k[1] == 0)
                {
                    continue;
                }

                //#Fucose <= #HexNAc. One Fucose modify one 
                if (k[4]!= 0 && k[4] > k[1] )
                {
                    continue;
                }

                //#NeuAc * 2 >= #Acetylation. One NeuAc can be modified with two Acetylation
                if (k[9]!= 0 && k[2]*2 < k[9])
                {
                    continue;
                }

                var ionMass = Glycan.GetMass(k);
                GlycanIon glycanIon = new GlycanIon(null, ionMass, k, glycan_mass - ionMass);
                glycanIons.Add(glycanIon);
            }

            return glycanIons.OrderBy(p=>p.IonMass).ToList();
        }

        /// <summary>
        /// Try to create all possible combinations from the glycan kind[]. And store the combination array in the _kinds list.
        /// </summary>
        /// <param name="kind"> ex. [2,2,0]</param>
        /// <param name="_kinds"></param>
        /// <param name="_keys"></param>
        private static void _GetCombinations(byte[] kind, List<byte[]> _kinds, HashSet<string> _keys) 
        {                                                                                            
            if (kind.Sum(p=>p) == 0)                                                                 
            {
                return; // if we don't have any monosaccharide, no need to generate the child ions.
            }
            else
            {
                for (int i = 0; i < kind.Length; i++) //traverse the kind array
                {
                    if (kind[i] >= 1)
                    {
                        byte[] akind = (byte[])kind.Clone();
                        akind[i]--;
                        if (akind.Sum(p => p) != 0)
                        {
                            if (!_keys.Contains(Glycan.GetKindString(akind)))
                            {
                                _keys.Add(Glycan.GetKindString(akind));
                                _kinds.Add((byte[])akind.Clone());
                                _GetCombinations(akind, _kinds, _keys);
                            }
                           
                        }

                    }
                }
            }
        }
    }
}
