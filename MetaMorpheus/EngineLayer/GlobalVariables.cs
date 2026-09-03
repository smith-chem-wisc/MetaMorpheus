global using obo = Omics.Modifications.IO.obo;
using Chemistry;
using Easy.Common.Extensions;
using EngineLayer.GlycoSearch;
using MassSpectrometry;
using Nett;
using Omics.Modifications;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using Omics.Modifications.IO;
using TopDownProteomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;
using System.Security.Cryptography;

namespace EngineLayer
{
    public static class GlobalVariables
    {
        public static string DecoyIdentifier { get; set; } = "DECOY";
        // for now, these are only used for error-checking in the command-line version.
        // compressed versions of the protein databases (e.g., .xml.gz) are also supported
        public static List<string> AcceptedDatabaseFormats { get; private set; }
        public static List<string> AcceptedSpectraFormats { get; private set; }

        private static List<Modification> _AllModsKnown;
        private static List<Modification> _AllRnaModsKnown;
        private static HashSet<string> _AllModTypesKnown;
        private static HashSet<string> _AllRnaModTypesKnown;
        private static List<Crosslinker> _KnownCrosslinkers;
        public static List<Modification> ProteaseMods = new List<Modification>();


        //Characters that aren't amino acids, but are reserved for special uses (motifs, delimiters, mods, etc)
        private static char[] _InvalidAminoAcids;

        // this affects output labels, etc. and can be changed to "Proteoform" for top-down searches
        public static AnalyteType AnalyteType;

        public static List<string> ErrorsReadingMods;

        /// <summary>
        /// Non-fatal things the user should know about that were noticed while starting up, and that are
        /// not about reading a modification. Surfaced beside <see cref="ErrorsReadingMods"/> by both front
        /// ends. Kept separate because that list's name is a promise about what is in it.
        /// </summary>
        public static List<string> StartupWarnings { get; private set; } = new List<string>();

        // mzLib keeps these as private constants, so the names are repeated rather than referenced.
        // They are the shipped files whose banner and header row seed the custom counterparts.
        private const string EmbeddedProteasesResourceName = "Proteomics.ProteolyticDigestion.proteases.tsv";
        private const string EmbeddedRnasesResourceName = "Transcriptomics.Digestion.rnases.tsv";

        /// <summary>Template seeded into the user's OGlycan_Custom.gdb. Embedded in EngineLayer.</summary>
        private const string EmbeddedCustomOGlycanResourceName = "EngineLayer.Glycan_Mods.OGlycan_Custom.gdb";

        /// <summary>Template seeded into the user's NGlycan_Custom.gdb. Embedded in EngineLayer.</summary>
        private const string EmbeddedCustomNGlycanResourceName = "EngineLayer.Glycan_Mods.NGlycan_Custom.gdb";

        /// <summary>Template seeded into Mods\CustomModifications.txt and Mods\RnaCustomModifications.txt.</summary>
        /// <remarks>
        /// The first line is the title CustomModWindow writes when it creates the file itself, so the GUI's
        /// append path and this template agree. The rest are '#' banner lines, which the modification format
        /// treats as comments -- see the shipped Mods.txt, which opens the same way.
        /// </remarks>
        private static string CustomModificationsTemplate(string analyte) =>
            "Custom Modifications" + Environment.NewLine +
            "################################## " + analyte + " modifications you add are stored here." + Environment.NewLine +
            "################################## Modifications added through the GUI are appended below this banner." + Environment.NewLine +
            "################################## One entry per modification, terminated by a line containing only //" + Environment.NewLine +
            "##################################   ID   <name>            the modification's name" + Environment.NewLine +
            "##################################   TG   <residues>        target, e.g. S or T or Y" + Environment.NewLine +
            "##################################   PP   <location>        Anywhere. / N-terminal. / C-terminal. / Peptide N-terminal. / Peptide C-terminal." + Environment.NewLine +
            "##################################   CF   <formula>         chemical formula, e.g. H1 N1 O2" + Environment.NewLine +
            "##################################   MM   <mass>            monoisotopic mass, if no formula is given" + Environment.NewLine +
            "##################################   MT   <type>            the group it is listed under" + Environment.NewLine +
            "################################## See Mods.txt in this folder for worked examples." + Environment.NewLine;

        // File locations
        public static string DataDir { get; private set; }
        public static string UserSpecifiedDataDir { get; set; }
        public static string CustomProteasePath => Path.Combine(DataDir, "proteases_custom.tsv");
        public static string CustomRnasePath => Path.Combine(DataDir, "rnase_custom.tsv");
        public static string CustomMonosaccharidePath => Path.Combine(DataDir, "MonosaccharidesCustom.tsv");

        /// <summary>
        /// The user's own O-glycan database, offered in the GlycoSearch task beside the shipped ones.
        /// </summary>
        /// <remarks>
        /// At the DataDir root rather than under Glycan_Mods\OGlycan\, for the same reason
        /// MonosaccharidesCustom.tsv is: Product.wxs gives Glycan_Mods and both of its subfolders a
        /// &lt;RemoveFolder On="both"/&gt;, so the installer owns those folders and a user's file in one of
        /// them is not somewhere we should be putting it. The cost is that it is not picked up by the
        /// directory sweep in LoadGlycans and has to be added to OGlycanDatabasePaths by name.
        /// </remarks>
        public static string CustomOGlycanDatabasePath => Path.Combine(DataDir, "OGlycan_Custom.gdb");

        /// <summary>
        /// The user's own N-glycan database, offered in the GlycoSearch task beside the shipped ones.
        /// At the DataDir root for the same reason as <see cref="CustomOGlycanDatabasePath"/>.
        /// </summary>
        public static string CustomNGlycanDatabasePath => Path.Combine(DataDir, "NGlycan_Custom.gdb");

        public static bool StopLoops { get; set; }
        public static string MetaMorpheusVersion { get; private set; }
        public static GlobalSettings GlobalSettings { get; set; }
        public static IEnumerable<Modification> UnimodDeserialized { get; private set; }
        public static IEnumerable<Modification> UniprotDeseralized { get; private set; }
        public static obo PsiModDeserialized { get; private set; }
        public static IEnumerable<Modification> AllModsKnown { get { return _AllModsKnown.AsEnumerable(); } }
        public static IEnumerable<Modification> AllRnaModsKnown { get { return _AllRnaModsKnown.AsEnumerable(); } }
        public static IEnumerable<string> AllModTypesKnown { get { return _AllModTypesKnown.AsEnumerable(); } }
        public static IEnumerable<string> AllRnaModTypesKnown { get { return _AllRnaModTypesKnown.AsEnumerable(); } }
        public static Dictionary<string, Modification> AllModsKnownDictionary { get; private set; }
        public static Dictionary<string, Modification> AllRnaModsKnownDictionary { get; private set; }
        public static Dictionary<string, string> AvailableUniProtProteomes { get; private set; }
        public static Dictionary<string, DissociationType> AllSupportedDissociationTypes { get; private set; }
        public static List<string> SeparationTypes { get; private set; }
        public static string ExperimentalDesignFileName { get; private set; }
        public static IEnumerable<Crosslinker> Crosslinkers { get { return _KnownCrosslinkers.AsEnumerable(); } }
        public static IEnumerable<char> InvalidAminoAcids { get { return _InvalidAminoAcids.AsEnumerable(); } }
        public static List<string> OGlycanDatabasePaths { get; private set; }
        public static List<string> NGlycanDatabasePaths { get; private set; }

        public static void SetUpGlobalVariables()
        {
            AcceptedDatabaseFormats = new List<string> { ".fasta", ".fa", ".xml", ".msp", ".msl" };
            // ".d" is the Bruker acquisition folder; the rest of the Bruker entries are the inner files a user may hand
            // us instead, which BrukerDataDirectory redirects to their parent ".d". Keep this list lower-case: every
            // consumer calls ToLowerInvariant() before Contains().
            AcceptedSpectraFormats = new List<string> { ".raw", ".mzml", ".mgf", ".msalign", ".baf", ".tdf", ".tdf_bin", ".tsf", ".tsf_bin", ".d" };
            AnalyteType = AnalyteType.Peptide;
            _InvalidAminoAcids = new char[] { 'X', 'B', 'J', 'Z', ':', '|', ';', '[', ']', '{', '}', '(', ')', '+', '-' };
            ExperimentalDesignFileName = "ExperimentalDesign.tsv";
            SeparationTypes = new List<string> { { "HPLC" }, { "CZE" } };

            StartupWarnings = new List<string>();

            SetMetaMorpheusVersion();
            SetUpDataDirectory();
            LoadCrosslinkers();
            LoadModifications();
            LoadRnaModifications();
            LoadGlycans();
            LoadCustomAminoAcids();
            SetUpGlobalSettings();
            LoadDissociationTypes();
            LoadAvailableProteomes();
            LoadDigestionAgents();
        }

        public static void AddMods(IEnumerable<Modification> modifications, bool modsAreFromTheTopOfProteinXml, bool isRna = false)
        {
            var allMods = isRna ? _AllRnaModsKnown : _AllModsKnown;
            var modTypes = isRna ? _AllRnaModTypesKnown : _AllModTypesKnown;

            foreach (var mod in modifications)
            {
                if (string.IsNullOrEmpty(mod.ModificationType) || string.IsNullOrEmpty(mod.IdWithMotif))
                {
                    ErrorsReadingMods.Add(mod.ToString() + Environment.NewLine + " has null or empty modification type");
                    continue;
                }
                if (allMods.Any(b => b.IdWithMotif.Equals(mod.IdWithMotif) && b.ModificationType.Equals(mod.ModificationType) && !b.Equals(mod)))
                {
                    if (modsAreFromTheTopOfProteinXml)
                    {
                        allMods.RemoveAll(p => p.IdWithMotif.Equals(mod.IdWithMotif) && p.ModificationType.Equals(mod.ModificationType) && !p.Equals(mod));
                        allMods.Add(mod);
                        modTypes.Add(mod.ModificationType);
                    }
                    else
                    {
                        ErrorsReadingMods.Add("Modification id and type are equal, but some fields are not! " +
                            "The following mod was not read in: " + Environment.NewLine + mod.ToString());
                    }
                    continue;
                }
                if (allMods.Any(b => b.IdWithMotif.Equals(mod.IdWithMotif) && b.ModificationType.Equals(mod.ModificationType)))
                {
                    // same ID, same mod type, and same mod properties; continue and don't output an error message
                    // this could result from reading in an XML database with mods annotated at the top
                    // that are already loaded in MetaMorpheus
                    continue;
                }
                if (allMods.Any(m => m.IdWithMotif == mod.IdWithMotif))
                {
                    // same ID but different mod types. This can happen if the user names a mod the same as a UniProt mod
                    // this is problematic because if a mod is annotated in the database, all we have to go on is an ID ("description" tag).
                    // so we don't know which mod to use, causing unnecessary ambiguity
                    if (modsAreFromTheTopOfProteinXml)
                    {
                        allMods.RemoveAll(p => p.IdWithMotif.Equals(mod.IdWithMotif) && !p.Equals(mod));
                        allMods.Add(mod);
                        modTypes.Add(mod.ModificationType);
                    }
                    else if (!mod.ModificationType.Equals("Unimod"))
                    {
                        ErrorsReadingMods.Add("Duplicate mod IDs! Skipping " + mod.ModificationType + ":" + mod.IdWithMotif);
                    }
                    continue;
                }

                // no errors! add the mod
                allMods.Add(mod);
                modTypes.Add(mod.ModificationType);
            }
        }

        public static void AddCrosslinkers(IEnumerable<Crosslinker> crosslinkers)
        {
            foreach (var linker in crosslinkers)
            {
                _KnownCrosslinkers.Add(linker);
            }
        }

        public static string CheckLengthOfOutput(string psmString)
        {
            if (psmString.Length > 32000 && GlobalSettings.WriteExcelCompatibleTSVs)
            {
                return "Output too long for Excel";
            }
            else
            {
                return psmString;
            }
        }

        public static void LoadCustomAminoAcids()
        {
            //read in all the amino acids (they already exist in mzlib, but there might be synthetic amino acids that need to be included)
            string aminoAcidPath = Path.Combine(DataDir, @"CustomAminoAcids", @"CustomAminoAcids.txt");
            if (File.Exists(aminoAcidPath)) //if it already exists
            {
                string[] aminoAcidLines = File.ReadAllLines(aminoAcidPath);
                List<Residue> residuesToAdd = new List<Residue>();
                for (int i = 1; i < aminoAcidLines.Length; i++)
                {

                    string[] line = aminoAcidLines[i].Split('\t').ToArray(); //tsv Name, one letter, monoisotopic, chemical formula
                    if (line.Length >= 4) //check something is there (not a blank line)
                    {
                        char letter = line[1][0];
                        if (InvalidAminoAcids.Contains(letter))
                        {
                            throw new MetaMorpheusException("Error while reading 'CustomAminoAcids.txt'. Line " + (i + 1).ToString() + " contains an invalid amino acid. (Ex: " + string.Join(", ", InvalidAminoAcids.Select(x => x.ToString())) + ")");
                        }
                        try
                        {
                            ChemicalFormula formula = ChemicalFormula.ParseFormula(line[3]);

                            //if it doesn't already exist or it does exist but has a different mass, add the entry
                            if (!(Residue.TryGetResidue(letter, out Residue residue))
                                || !(formula.Formula.Equals(residue.ThisChemicalFormula.Formula)))
                            {
                                residuesToAdd.Add(new Residue(line[0], letter, line[1], formula, ModificationSites.Any));
                            }
                        }
                        catch
                        {
                            throw new MetaMorpheusException("Error while reading 'CustomAminoAcids.txt'. Line " + (i + 1).ToString() + " was not in the correct format.");
                        }
                    }
                }
                Residue.AddNewResiduesToDictionary(residuesToAdd);
            }
            else //create it so that it can be manipulated
            {
                WriteAminoAcidsFile();
            }
        }

        public static void WriteAminoAcidsFile()
        {
            string directory = Path.Combine(DataDir, @"CustomAminoAcids");
            if (!Directory.Exists(directory))
            {
                Directory.CreateDirectory(directory);
            }
            string aminoAcidPath = Path.Combine(DataDir, @"CustomAminoAcids", @"CustomAminoAcids.txt");
            List<string> linesToWrite = new List<string> { "Name\tOneLetterAbbr.\tMonoisotopicMass\tChemicalFormula" };
            for (char letter = 'A'; letter <= 'Z'; letter++) //just the basic residues
            {
                if (Residue.TryGetResidue(letter, out Residue residue))
                {
                    linesToWrite.Add(residue.Name + '\t' + residue.Letter.ToString() + '\t' + residue.MonoisotopicMass.ToString() + '\t' + residue.ThisChemicalFormula.Formula);
                }
            }
            File.WriteAllLines(aminoAcidPath, linesToWrite.ToArray());
        }

        // Does the same thing as Process.Start() except it works on .NET Core
        public static void StartProcess(string path, bool useNotepadToOpenToml = false)
        {
            var p = new Process();

            p.StartInfo = new ProcessStartInfo()
            {
                UseShellExecute = true,
                FileName = path
            };

            if (useNotepadToOpenToml && Path.GetExtension(path).ToLowerInvariant() == ".toml" && RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
            {
                p.StartInfo.FileName = "notepad.exe";
                p.StartInfo.Arguments = path;
            }

            p.Start();
        }

        public static void CopyFilesRecursively(DirectoryInfo source, DirectoryInfo target)
        {
            //https://stackoverflow.com/questions/58744/copy-the-entire-contents-of-a-directory-in-c-sharp
            foreach (DirectoryInfo dir in source.GetDirectories())
            {
                CopyFilesRecursively(dir, target.CreateSubdirectory(dir.Name));
            }
            foreach (FileInfo file in source.GetFiles())
            {
                file.CopyTo(Path.Combine(target.FullName, file.Name));
            }
        }

        /// <summary>
        /// Gets the file extension, with the option to keep .gz appended for compressed files
        /// </summary>
        public static string GetFileExtension(string fileWithExtension, bool getUncompressedExtension = true)
        {
            string extension = string.Empty;
            StringBuilder sb = new StringBuilder();

            for (int i = fileWithExtension.Length - 1; i >= 0; i--)
            {
                char c = fileWithExtension[i];

                sb.Append(c);

                if (c == '.')
                {
                    extension = new string(sb.ToString().Reverse().ToArray());

                    if (!extension.ToLowerInvariant().EndsWith("gz") || extension.Count(p => p == '.') >= 2)
                    {
                        break;
                    }
                }
            }

            if (getUncompressedExtension && extension.ToLowerInvariant().EndsWith("gz"))
            {
                int indexOfGz = extension.ToLowerInvariant().IndexOf("gz");

                for (int i = indexOfGz; i >= 0; i--)
                {
                    if (extension[i] == '.')
                    {
                        extension = extension.Substring(0, i);
                        break;
                    }
                }
            }

            return extension;
        }

        public static string GetFilenameWithoutExtension(string path)
        {
            Path.GetFileNameWithoutExtension("");
            var filename = Path.GetFileName(path);
            string extension = GetFileExtension(filename, getUncompressedExtension: false);

            if (extension == string.Empty)
            {
                return filename;
            }

            return filename.Replace(extension, string.Empty);
        }

        private static void SetMetaMorpheusVersion()
        {
            // get version of this MetaMorpheus instance
            MetaMorpheusVersion = typeof(GlobalVariables).Assembly.GetName().Version.ToString();

            if (MetaMorpheusVersion.Equals("1.0.0.0"))
            {
#if DEBUG
                MetaMorpheusVersion = "Not a release version. DEBUG.";
#else
                MetaMorpheusVersion = "Not a release version.";
#endif
            }
            else
            {
                // as of 0.0.277, AppVeyor appends the build number
                // this is intentional; it's to avoid conflicting AppVeyor build numbers
                // trim the build number off the version number for displaying/checking versions, etc
                var foundIndexes = new List<int>();
                for (int i = 0; i < MetaMorpheusVersion.Length; i++)
                {
                    if (MetaMorpheusVersion[i] == '.')
                    {
                        foundIndexes.Add(i);
                    }
                }
                MetaMorpheusVersion = MetaMorpheusVersion.Substring(0, foundIndexes.Last());
            }
        }

        private static void SetUpDataDirectory()
        {
            // get data directory
            var pathToProgramFiles = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
            if (!String.IsNullOrWhiteSpace(pathToProgramFiles) && AppDomain.CurrentDomain.BaseDirectory.Contains(pathToProgramFiles)
                && !AppDomain.CurrentDomain.BaseDirectory.Contains("Jenkins"))
            {
                DataDir = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "MetaMorpheus");
            }
            else
            {
                DataDir = AppDomain.CurrentDomain.BaseDirectory;
            }

            if (UserSpecifiedDataDir != null)
            {
                if (!Directory.Exists(UserSpecifiedDataDir))
                {
                    CopyFilesRecursively(new DirectoryInfo(DataDir), new DirectoryInfo(UserSpecifiedDataDir));
                }

                DataDir = UserSpecifiedDataDir;
            }
        }

        private static void LoadCrosslinkers()
        {
            _KnownCrosslinkers = new List<Crosslinker>();

            // load default crosslinkers
            string crosslinkerLocation = Path.Combine(DataDir, @"Data", @"Crosslinkers.tsv");
            AddCrosslinkers(Crosslinker.LoadCrosslinkers(crosslinkerLocation));

            // load custom crosslinkers
            string customCrosslinkerLocation = Path.Combine(DataDir, @"Data", @"CustomCrosslinkers.tsv");

            // Header row only, with no banner: LoadCrosslinkers skips line 1 and parses every line after
            // it, so a comment banner here would be read as a crosslinker and throw on its columns.
            CustomDataFile.EnsureExists(customCrosslinkerLocation,
                () => CustomDataFile.BannerAndHeaderFromFile(crosslinkerLocation, "Name\t"),
                "custom crosslinker");

            if (File.Exists(customCrosslinkerLocation))
            {
                AddCrosslinkers(Crosslinker.LoadCrosslinkers(customCrosslinkerLocation));
            }
        }

        private static void LoadModifications()
        {
            _AllModsKnown = new List<Modification>();
            _AllModTypesKnown = new HashSet<string>();
            ErrorsReadingMods = new List<string>();
            AllModsKnownDictionary = new Dictionary<string, Modification>();

            UnimodDeserialized = Loaders.LoadUnimod(Path.Combine(DataDir, @"Data", @"unimod.xml")).ToList();
            PsiModDeserialized = Loaders.LoadPsiMod(Path.Combine(DataDir, @"Data", @"PSI-MOD.obo.xml"));
            var formalChargesDictionary = Loaders.GetFormalChargesDictionary(PsiModDeserialized);
            UniprotDeseralized = Loaders.LoadUniprot(Path.Combine(DataDir, @"Data", @"ptmlist.txt"), formalChargesDictionary).ToList();

            // Seeded before the sweep below picks it up. The template is a title line plus a '#' banner,
            // so it contributes no modifications until the user or the GUI adds one.
            CustomDataFile.EnsureExists(Path.Combine(DataDir, @"Mods", "CustomModifications.txt"),
                () => CustomModificationsTemplate("Protein"), "custom modification");

            foreach (var modFile in Directory.GetFiles(Path.Combine(DataDir, @"Mods")))
            {
                if (modFile.Contains("glyco.txt"))
                {
                    // Glycan modifications are handled separately in LoadGlycans()
                    continue;
                }
                if (modFile.Contains("Rna"))
                    continue;
                AddMods(ModificationLoader.ReadModsFromFile(modFile, out var errorMods), false);
            }

            AddMods(UniprotDeseralized.OfType<Modification>(), false);
            AddMods(UnimodDeserialized.OfType<Modification>(), false);
            
            foreach (Modification mod in AllModsKnown)
            {
                if (!AllModsKnownDictionary.ContainsKey(mod.IdWithMotif))
                {
                    AllModsKnownDictionary.Add(mod.IdWithMotif, mod);
                }
                // no error thrown if multiple mods with this ID are present - just pick one
            }
        }

        private static void LoadRnaModifications()
        {
            _AllRnaModsKnown = new List<Modification>();
            _AllRnaModTypesKnown = new HashSet<string>();
            AllRnaModsKnownDictionary = new Dictionary<string, Modification>();

            // RNA Mods is an embedded resources: It gets packed into the DLL so we do not need to worry about the installer. 
            var assembly = typeof(GlobalVariables).Assembly;
            var resourceName = "EngineLayer.Mods.RnaMods.txt";

            using (var stream = assembly.GetManifestResourceStream(resourceName))
            using (var reader = new StreamReader(stream))
            {
                string fileContent = reader.ReadToEnd();
                var mods = ModificationLoader.ReadModsFromString(fileContent, out var errors);
                AddMods(mods, false, true);
            }

            var customModsPath = Path.Combine(DataDir, @"Mods", "RnaCustomModifications.txt");
            CustomDataFile.EnsureExists(customModsPath,
                () => CustomModificationsTemplate("RNA"), "custom RNA modification");
            if (File.Exists(customModsPath))
            {
                AddMods(ModificationLoader.ReadModsFromFile(customModsPath, out var errorMods), false, true);
            }

            // populate mod types and dictionary
            _AllRnaModsKnown.Select(mod => mod.ModificationType)
                .Distinct()
                .ForEach(type => _AllRnaModTypesKnown.Add(type));

            AllRnaModsKnownDictionary = _AllRnaModsKnown.Where(p => p.OriginalId != "").ToDictionary(p => $"{p.IdWithMotif}");
        }

        private static void LoadGlycans()
        {
            // Custom monosaccharides must be registered FIRST so any custom tokens are recognized
            // by the glycan-database parsers below. EnsureCustomMonosaccharideFileExists seeds the
            // file (from the embedded template, or a carried-over legacy copy) if it's missing, so
            // LoadCustomMonosaccharides always has a file to read here.
            GlycanDatabase.EnsureCustomMonosaccharideFileExists(CustomMonosaccharidePath);
            GlycanDatabase.LoadCustomMonosaccharides(CustomMonosaccharidePath);

            // Seed the user's own database before anything reads it. It is header-less and its template is
            // all comment lines, so a freshly seeded file contributes no glycans and the two steps below are
            // independent of it. See CustomDataFile for the recipe every custom file follows.
            CustomDataFile.EnsureExists(CustomOGlycanDatabasePath,
                () => CustomDataFile.EmbeddedText(typeof(GlobalVariables).Assembly, EmbeddedCustomOGlycanResourceName),
                "custom O-glycan database");

            CustomDataFile.EnsureExists(CustomNGlycanDatabasePath,
                () => CustomDataFile.EmbeddedText(typeof(GlobalVariables).Assembly, EmbeddedCustomNGlycanResourceName),
                "custom N-glycan database");

            OGlycanDatabasePaths = new List<string>();
            NGlycanDatabasePaths = new List<string>();

            foreach (var glycanFile in Directory.GetFiles(Path.Combine(DataDir, @"Glycan_Mods", @"OGlycan")))
            {
                OGlycanDatabasePaths.Add(glycanFile);
            }

            // Added by name because it deliberately does not live in the swept folder -- see
            // CustomOGlycanDatabasePath. Guarded so that a seeding failure earlier cannot put a path to a
            // file that is not there into the list the task window offers.
            if (File.Exists(CustomOGlycanDatabasePath))
            {
                OGlycanDatabasePaths.Add(CustomOGlycanDatabasePath);
            }

            foreach (var glycanFile in Directory.GetFiles(Path.Combine(DataDir, @"Glycan_Mods", @"NGlycan")))
            {
                NGlycanDatabasePaths.Add(glycanFile);
            }

            if (File.Exists(CustomNGlycanDatabasePath))
            {
                NGlycanDatabasePaths.Add(CustomNGlycanDatabasePath);
            }

            //Add Glycan mod into AllModsKnownDictionary, currently this is for MetaDraw.
            //The reason why not include Glycan into modification database is for users to apply their own database.
            foreach (var path in OGlycanDatabasePaths)
            {
                var oGlycans = GlycanDatabase.LoadGlycan(path, false, true);
                foreach (var glycan in oGlycans)
                {
                    if (!AllModsKnownDictionary.ContainsKey(glycan.IdWithMotif))
                    {
                        AllModsKnownDictionary.Add(glycan.IdWithMotif, glycan);
                    }
                    _AllModsKnown.Add(glycan);
                }
            }
            foreach (var path in NGlycanDatabasePaths)
            {
                var nGlycans = GlycanDatabase.LoadGlycan(path, false, false);
                foreach (var glycan in nGlycans)
                {
                    if (!AllModsKnownDictionary.ContainsKey(glycan.IdWithMotif))
                    {
                        AllModsKnownDictionary.Add(glycan.IdWithMotif, glycan);
                    }
                    _AllModsKnown.Add(glycan);
                }
            }
            LoadTxtGlycan();
        }

        private static void LoadDissociationTypes()
        {
            // set up dissociation types
            AllSupportedDissociationTypes = new Dictionary<string, DissociationType> {
                { DissociationType.CID.ToString(), DissociationType.CID },
                { DissociationType.ECD.ToString(), DissociationType.ECD },
                { DissociationType.ETD.ToString(), DissociationType.ETD },
                { DissociationType.HCD.ToString(), DissociationType.HCD },
                { DissociationType.EThcD.ToString(), DissociationType.EThcD },
                { DissociationType.Custom.ToString(), DissociationType.Custom },
                { DissociationType.LowCID.ToString(), DissociationType.LowCID},

                // allow reading from scan header (autodetect dissociation type)
                { DissociationType.Autodetect.ToString(), DissociationType.Autodetect}
            };
        }

        private static void LoadAvailableProteomes()
        {
            AvailableUniProtProteomes = ProteinDbRetriever.UniprotProteomesList(Path.Combine(DataDir,@"Proteomes",@"availableUniProtProteomes.txt.gz"));
        }
        private static void SetUpGlobalSettings()
        {
            // save/load settings
            string settingsPath = Path.Combine(DataDir, @"settings.toml");
            if (!File.Exists(settingsPath) && !new DirectoryInfo(DataDir).Attributes.HasFlag(FileAttributes.ReadOnly))
            {
                Toml.WriteFile<GlobalSettings>(new GlobalSettings(), settingsPath);
            }

            if (File.Exists(settingsPath))
            {
                GlobalSettings = Toml.ReadFile<GlobalSettings>(settingsPath);
            }
        }

        /// <summary>
        /// Convert glyco.txt into Glycan objects and add them to AllModsKnown.
        /// </summary>
        private static void LoadTxtGlycan()
        {
            string glycoFile = Path.Combine(DataDir, @"Mods", "glyco.txt");
            var glycoMods = ModificationLoader.ReadModsFromFile(glycoFile, out var errorMods);
            foreach (var glycoMod in glycoMods)
            {
                var kind = GlycanDatabase.String2Kind(glycoMod.OriginalId);

                // If we cannot parse the glycan string, we add the glycoMod as a normal modification.
                if (kind.Sum(p => p) == 0)
                {
                    _AllModsKnown.Add(glycoMod);
                    continue;
                }

                Glycan glycan;
                if (glycoMod.ModificationType == "N-linked glycosylation")
                {
                    glycan = new Glycan(kind, glycoMod.Target.ToString(), GlycanType.N_glycan);
                    glycan.Ions = GlycanDatabase.OGlycanCompositionCombinationChildIons(kind);
                }
                else
                {
                    glycan = new Glycan(kind, glycoMod.Target.ToString(), GlycanType.O_glycan);
                    glycan.Ions = GlycanDatabase.OGlycanCompositionCombinationChildIons(kind);
                }
                _AllModsKnown.Add(glycan);
            }
        }

        private static void LoadDigestionAgents()
        {
            // Seed first, then load: the template is header-only, so a freshly seeded file contributes
            // nothing and the two steps are independent. See CustomDataFile for the recipe every custom
            // file follows.
            CustomDataFile.EnsureExists(CustomProteasePath,
                () => CustomDataFile.BannerAndHeaderFrom(typeof(ProteaseDictionary).Assembly,
                    EmbeddedProteasesResourceName, "Name\t"),
                "custom protease");

            CustomDataFile.EnsureExists(CustomRnasePath,
                () => CustomDataFile.BannerAndHeaderFrom(typeof(RnaseDictionary).Assembly,
                    EmbeddedRnasesResourceName, "Name\t"),
                "custom rnase");

            if (File.Exists(CustomProteasePath))
            {
                try
                {
                    var mods = ProteaseDictionary.LoadEmbeddedProteaseMods();
                    var result = ProteaseDictionary.LoadAndMergeCustomProteases(CustomProteasePath, mods);
                    ReportSkippedCustomEntries(result.Skipped, "protease", CustomProteasePath);
                }
                catch (Exception e)
                {
                    throw new MetaMorpheusException($"Error loading custom proteases with error message: {e.Message}", e);
                }
            }

            if (File.Exists(CustomRnasePath))
            {
                try
                {
                    var result = RnaseDictionary.LoadAndMergeCustomRnases(CustomRnasePath);
                    ReportSkippedCustomEntries(result.Skipped, "rnase", CustomRnasePath);
                }
                catch (Exception e)
                {
                    throw new MetaMorpheusException($"Error loading custom rnases with error message: {e.Message}", e);
                }
            }
        }

        /// <summary>
        /// mzLib refuses to let a custom digestion agent shadow one of its own, and reports the collision
        /// through <c>CustomDigestionAgentLoadResult.Skipped</c> rather than throwing, specifically so the
        /// caller can tell the user. Nothing consumed that before, so a user who named a custom protease
        /// "trypsin" got silence and a protease that was not theirs.
        /// </summary>
        private static void ReportSkippedCustomEntries(IReadOnlyList<string> skipped, string kind, string path)
        {
            if (skipped == null || skipped.Count == 0)
            {
                return;
            }

            StartupWarnings.Add($"{skipped.Count} custom {kind}(s) in {Path.GetFileName(path)} were ignored because "
                + $"a built-in {kind} already uses the same name: {string.Join(", ", skipped.Select(p => "'" + p + "'"))}. "
                + $"Rename them in {path} if you meant to define your own.");
        }
    }
}
