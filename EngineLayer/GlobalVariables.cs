using Chemistry;
using MassSpectrometry;
using Nett;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using UsefulProteomicsDatabases;

namespace EngineLayer
{
    public static class GlobalVariables
    {
        // for now, these are only used for error-checking in the command-line version.
        // compressed versions of the protein databases (e.g., .xml.gz) are also supported
        public static List<string> AcceptedDatabaseFormats { get; private set; }
        public static List<string> AcceptedSpectraFormats { get; private set; }

        private static List<Modification> _AllModsKnown;
        private static HashSet<string> _AllModTypesKnown;
        private static List<Crosslinker> _KnownCrosslinkers;
        public static List<Modification> ProteaseMods = new List<Modification>();


        //Characters that aren't amino acids, but are reserved for special uses (motifs, delimiters, mods, etc)
        private static char[] _InvalidAminoAcids;

        // this affects output labels, etc. and can be changed to "Proteoform" for top-down searches
        public static string AnalyteType;

        public static List<string> ErrorsReadingMods;

        // File locations
        public static string DataDir { get; private set; }
        public static string UserSpecifiedDataDir { get; set; }

        public static bool StopLoops { get; set; }
        public static string MetaMorpheusVersion { get; private set; }
        public static GlobalSettings GlobalSettings { get; set; }
        public static IEnumerable<Modification> UnimodDeserialized { get; private set; }
        public static IEnumerable<Modification> UniprotDeseralized { get; private set; }
        public static UsefulProteomicsDatabases.Generated.obo PsiModDeserialized { get; private set; }
        public static IEnumerable<Modification> AllModsKnown { get { return _AllModsKnown.AsEnumerable(); } }
        public static IEnumerable<string> AllModTypesKnown { get { return _AllModTypesKnown.AsEnumerable(); } }
        public static Dictionary<string, Modification> AllModsKnownDictionary { get; private set; }
        public static Dictionary<string, string> AvailableUniProtProteomes { get; private set; }
        public static Dictionary<string, DissociationType> AllSupportedDissociationTypes { get; private set; }
        public static List<string> SeparationTypes { get; private set; }
        public static string ExperimentalDesignFileName { get; private set; }
        public static IEnumerable<Crosslinker> Crosslinkers { get { return _KnownCrosslinkers.AsEnumerable(); } }
        public static IEnumerable<char> InvalidAminoAcids { get { return _InvalidAminoAcids.AsEnumerable(); } }
        public static List<string> OGlycanLocations { get; private set; }
        public static List<string> NGlycanLocations { get; private set; }

        public static void SetUpGlobalVariables()
        {
            Loaders.LoadElements();
            AcceptedDatabaseFormats = new List<string> { ".fasta", ".fa", ".xml", ".msp" };
            AcceptedSpectraFormats = new List<string> { ".raw", ".mzml", ".mgf" };
            AnalyteType = "Peptide";
            _InvalidAminoAcids = new char[] { 'X', 'B', 'J', 'Z', ':', '|', ';', '[', ']', '{', '}', '(', ')', '+', '-' };
            ExperimentalDesignFileName = "ExperimentalDesign.tsv";
            SeparationTypes = new List<string> { { "HPLC" }, { "CZE" } };

            SetMetaMorpheusVersion();
            SetUpDataDirectory();
            LoadCrosslinkers();
            LoadModifications();
            LoadGlycans();
            LoadCustomAminoAcids();
            SetUpGlobalSettings();
            LoadDissociationTypes();
            LoadAvailableProteomes();
        }

        public static void AddMods(IEnumerable<Modification> modifications, bool modsAreFromTheTopOfProteinXml)
        {
            foreach (var mod in modifications)
            {
                if (string.IsNullOrEmpty(mod.ModificationType) || string.IsNullOrEmpty(mod.IdWithMotif))
                {
                    ErrorsReadingMods.Add(mod.ToString() + Environment.NewLine + " has null or empty modification type");
                    continue;
                }
                if (AllModsKnown.Any(b => b.IdWithMotif.Equals(mod.IdWithMotif) && b.ModificationType.Equals(mod.ModificationType) && !b.Equals(mod)))
                {
                    if (modsAreFromTheTopOfProteinXml)
                    {
                        _AllModsKnown.RemoveAll(p => p.IdWithMotif.Equals(mod.IdWithMotif) && p.ModificationType.Equals(mod.ModificationType) && !p.Equals(mod));
                        _AllModsKnown.Add(mod);
                        _AllModTypesKnown.Add(mod.ModificationType);
                    }
                    else
                    {
                        ErrorsReadingMods.Add("Modification id and type are equal, but some fields are not! " +
                            "The following mod was not read in: " + Environment.NewLine + mod.ToString());
                    }
                    continue;
                }
                else if (AllModsKnown.Any(b => b.IdWithMotif.Equals(mod.IdWithMotif) && b.ModificationType.Equals(mod.ModificationType)))
                {
                    // same ID, same mod type, and same mod properties; continue and don't output an error message
                    // this could result from reading in an XML database with mods annotated at the top
                    // that are already loaded in MetaMorpheus
                    continue;
                }
                else if (AllModsKnown.Any(m => m.IdWithMotif == mod.IdWithMotif))
                {
                    // same ID but different mod types. This can happen if the user names a mod the same as a UniProt mod
                    // this is problematic because if a mod is annotated in the database, all we have to go on is an ID ("description" tag).
                    // so we don't know which mod to use, causing unnecessary ambiguity
                    if (modsAreFromTheTopOfProteinXml)
                    {
                        _AllModsKnown.RemoveAll(p => p.IdWithMotif.Equals(mod.IdWithMotif) && !p.Equals(mod));
                        _AllModsKnown.Add(mod);
                        _AllModTypesKnown.Add(mod.ModificationType);
                    }
                    else if (!mod.ModificationType.Equals("Unimod"))
                    {
                        ErrorsReadingMods.Add("Duplicate mod IDs! Skipping " + mod.ModificationType + ":" + mod.IdWithMotif);
                    }
                    continue;
                }
                else
                {
                    // no errors! add the mod
                    _AllModsKnown.Add(mod);
                    _AllModTypesKnown.Add(mod.ModificationType);
                }
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

            foreach (var modFile in Directory.GetFiles(Path.Combine(DataDir, @"Mods")))
            {
                AddMods(PtmListLoader.ReadModsFromFile(modFile, out var errorMods), false);
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
            ProteaseMods = UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(Path.Combine(DataDir, @"Mods", @"ProteaseMods.txt"), out var errors).ToList();
            ProteaseDictionary.Dictionary = ProteaseDictionary.LoadProteaseDictionary(Path.Combine(DataDir, @"ProteolyticDigestion", @"proteases.tsv"), ProteaseMods);
        }

        private static void LoadGlycans()
        {
            OGlycanLocations = new List<string>();
            NGlycanLocations = new List<string>();

            foreach (var glycanFile in Directory.GetFiles(Path.Combine(DataDir, @"Glycan_Mods", @"OGlycan")))
            {
                OGlycanLocations.Add(glycanFile);
            }

            foreach (var glycanFile in Directory.GetFiles(Path.Combine(DataDir, @"Glycan_Mods", @"NGlycan")))
            {
                NGlycanLocations.Add(glycanFile);
            }

            //Add Glycan mod into AllModsKnownDictionary, currently this is for MetaDraw.
            //The reason why not include Glycan into modification database is for users to apply their own database.
            foreach (var path in OGlycanLocations)
            {
                var og = GlycanDatabase.LoadGlycan(path, false, false);
                foreach (var g in og)
                {
                    var ogmod = Glycan.OGlycanToModification(g);
                    if (!AllModsKnownDictionary.ContainsKey(ogmod.IdWithMotif))
                    {
                        AllModsKnownDictionary.Add(ogmod.IdWithMotif, ogmod);
                    }
                }
            }
            foreach (var path in NGlycanLocations)
            {
                var og = GlycanDatabase.LoadGlycan(path, false, false);
                foreach (var g in og)
                {
                    var ogmod = Glycan.OGlycanToModification(g);
                    if (!AllModsKnownDictionary.ContainsKey(ogmod.IdWithMotif))
                    {
                        AllModsKnownDictionary.Add(ogmod.IdWithMotif, ogmod);
                    }
                }
            }
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
    }
}