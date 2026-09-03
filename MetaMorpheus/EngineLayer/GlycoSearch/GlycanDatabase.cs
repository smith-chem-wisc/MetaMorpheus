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
                    if (!line.Contains("HexNAc"))  // use the first data line to determine the format (kind / structure) of glycan database.
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
        /// content, and the specific problem (collision with built-in, bad char, non-numeric mass).
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
                }
            }
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

            if (ContainsEntry(databasePath, entry))
            {
                throw new MetaMorpheusException($"Could not add the glycan: '{fileName}' already contains \"{entry}\".");
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
                ValidateStructureEntry(entry, fileName, isOGlycan);
                return GlycanLineFormat.Structure;
            }

            ValidateCompositionEntry(entry);
            return GlycanLineFormat.Composition;
        }

        private static void ValidateStructureEntry(string entry, string fileName, bool isOGlycan)
        {
            // The same character check the loader runs, so an unknown monosaccharide code is rejected here
            // rather than being silently counted as zero mass at search time.
            ValidateStructureLine(entry, fileName, 1);

            int depth = 0;
            foreach (char c in entry)
            {
                if (c == '(')
                {
                    depth++;
                }
                else if (c == ')')
                {
                    depth--;
                    if (depth < 0)
                    {
                        throw new MetaMorpheusException(
                            $"Could not add the glycan \"{entry}\": the parentheses do not balance -- a ')' closes a branch that was never opened.");
                    }
                }
            }
            if (depth != 0)
            {
                throw new MetaMorpheusException(
                    $"Could not add the glycan \"{entry}\": the parentheses do not balance -- {depth} branch(es) are left open.");
            }

            try
            {
                // Parsed exactly as LoadStructureGlycan would parse it: Struct2Glycan is the only thing that
                // knows whether the nesting describes a tree it can actually build.
                List<Glycan> parsed = Glycan.Struct2Glycan(entry, 1, isOGlycan);
                if (parsed == null || parsed.Count == 0)
                {
                    throw new MetaMorpheusException($"Could not add the glycan \"{entry}\": it did not parse into any glycan.");
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

        private static void ValidateCompositionEntry(string entry)
        {
            // Name(count), repeated -- e.g. HexNAc(2)Hex(5)NeuAc(1). Checked here rather than by handing it
            // to String2Kind, which throws KeyNotFoundException on an unknown name and quietly accepts a
            // trailing unclosed group.
            Match match = Regex.Match(entry, @"^(?:(?<name>[A-Za-z][A-Za-z0-9]*)\((?<count>\d+)\))+$");
            if (!match.Success)
            {
                throw new MetaMorpheusException(
                    $"Could not add the glycan \"{entry}\": it is neither a structure -- which starts with '(', e.g. (N(H(A))) -- " +
                    "nor a composition, which is a name and a count repeated, e.g. HexNAc(2)Hex(5).");
            }

            CaptureCollection names = match.Groups["name"].Captures;
            CaptureCollection counts = match.Groups["count"].Captures;
            HashSet<string> seen = new HashSet<string>(StringComparer.Ordinal);

            for (int i = 0; i < names.Count; i++)
            {
                string name = names[i].Value;
                if (!Glycan.NameCharDic.ContainsKey(name))
                {
                    string allowed = string.Join(", ", Glycan.NameCharDic.Keys);
                    throw new MetaMorpheusException(
                        $"Could not add the glycan \"{entry}\": '{name}' is not a monosaccharide MetaMorpheus knows. Known: {allowed}. " +
                        "A monosaccharide it does not ship with must be declared in MonosaccharidesCustom.tsv first.");
                }
                if (!seen.Add(name))
                {
                    throw new MetaMorpheusException(
                        $"Could not add the glycan \"{entry}\": '{name}' appears more than once. Give each monosaccharide a single total count.");
                }
                if (!byte.TryParse(counts[i].Value, NumberStyles.Integer, CultureInfo.InvariantCulture, out byte count))
                {
                    throw new MetaMorpheusException(
                        $"Could not add the glycan \"{entry}\": the count for '{name}' must be a whole number between 0 and 255.");
                }
            }
        }

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

        /// <summary>Whether the database already holds this exact glycan, so it is not added twice.</summary>
        private static bool ContainsEntry(string databasePath, string entry)
        {
            if (!File.Exists(databasePath))
            {
                return false;
            }

            foreach (string line in File.ReadLines(databasePath))
            {
                // Split on tab first: the shipped .txt databases carry name and mass columns after the
                // glycan, and it is the glycan that has to be unique.
                if (!IsCommentOrBlank(line) && line.Split('\t')[0].Trim().Equals(entry, StringComparison.Ordinal))
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Ensure the MonosaccharidesCustom.tsv exists in the directory. If the file is missing, 
        /// write the embedded full 85-line documented template—instructions, column spec, the built-in name/code table,
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
                while (lines.Peek() != -1)
                {
                    string rawLine = lines.ReadLine();

                    // Skipped explicitly rather than left to the Hex test below: a comment that quotes a
                    // composition -- which the documented template does, to show the format -- would otherwise
                    // reach String2Kind and die on a dictionary lookup naming neither the file nor the line.
                    if (IsCommentOrBlank(rawLine))
                    {
                        continue;
                    }

                    string line = rawLine.Split('\t').First();

                    if (!(line.Contains("HexNAc") || line.Contains("Hex"))) // Make sure the line is a glycan line. The line should contain HexNAc or Hex.
                    {
                        continue;
                    }

                    var kind = String2Kind(line);  // Convert the database string to kind[] format (byte array).

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
                    // and any database a user has annotated -- reaches ValidateStructureLine and throws on
                    // the '#' itself during startup, before any window opens.
                    if (IsCommentOrBlank(line))
                    {
                        continue;
                    }

                    // ValidateStructureLine catches characters the parser would otherwise silently
                    // miscount as zero-mass. It consults the live monosaccharide registry, so any
                    // codes registered via MonosaccharidesCustom.tsv are accepted here too.
                    ValidateStructureLine(line.Trim(), filePath, lineNumber);

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

        // Valid characters in structure-format lines: parens plus any single-char monosaccharide
        // code currently registered with Glycan (built-ins + customs loaded from
        // MonosaccharidesCustom.tsv). Struct2Glycan silently miscounts unknown chars (no entry in
        // CharMassDic), so we pre-validate to give the user a clear error instead of a silently
        // wrong glycan mass.
        private static void ValidateStructureLine(string trimmedLine, string filePath, int lineNumber)
        {
            foreach (char c in trimmedLine)
            {
                if (c == '(' || c == ')') continue;
                if (!Glycan.CharMassDic.ContainsKey(c))
                {
                    string allowed = string.Concat(Glycan.CharMassDic.Keys);
                    throw new MetaMorpheusException(
                        $"Could not parse glycan structure in '{Path.GetFileName(filePath)}' at line {lineNumber}: \"{trimmedLine}\". " +
                        $"Unrecognized character '{c}'. Allowed: parentheses and one of {allowed}. " +
                        "A monosaccharide MetaMorpheus does not ship with must be declared in MonosaccharidesCustom.tsv first.");
                }
            }
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
