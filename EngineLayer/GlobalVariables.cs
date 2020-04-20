using Chemistry;
using MassSpectrometry;
using Nett;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace EngineLayer
{
    public static class GlobalVariables
    {
        // for now, these are only used for error-checking in the command-line version.
        // compressed versions of the protein databases (e.g., .xml.gz) are also supported
        public static List<string> AcceptedDatabaseFormats = new List<string> { ".fasta", ".fa", ".xml" };
        public static List<string> AcceptedSpectraFormats = new List<string> { ".raw", ".mzml", ".mgf" };

        private static List<Modification> _AllModsKnown = new List<Modification>();
        private static HashSet<string> _AllModTypesKnown = new HashSet<string>();
        private static List<Crosslinker> _KnownCrosslinkers = new List<Crosslinker>();
        private static List<string> _SeparationTypes = new List<string>();

        //Characters that aren't amino acids, but are reserved for special uses (motifs, delimiters, mods, etc)
        private static char[] _InvalidAminoAcids = new char[] { 'X', 'B', 'J', 'Z', ':', '|', ';', '[', ']', '{', '}', '(', ')', '+', '-' };

        // this affects output labels, etc. and can be changed to "Proteoform" for top-down searches
        public static string AnalyteType = "Peptide"; 

        static GlobalVariables()
        {
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

            {
                var pathToProgramFiles = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
                if (!String.IsNullOrWhiteSpace(pathToProgramFiles) && AppDomain.CurrentDomain.BaseDirectory.Contains(pathToProgramFiles) && !AppDomain.CurrentDomain.BaseDirectory.Contains("Jenkins"))
                {
                    DataDir = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "MetaMorpheus");
                }
                else
                {
                    DataDir = AppDomain.CurrentDomain.BaseDirectory;
                }
            }

            ElementsLocation = Path.Combine(DataDir, @"Data", @"elements.dat");
            UsefulProteomicsDatabases.Loaders.LoadElements();

            AddSeparationTypes(new List<string> { { "HPLC" }, { "CZE" } });

            // load default crosslinkers
            string crosslinkerLocation = Path.Combine(DataDir, @"Data", @"Crosslinkers.tsv");
            AddCrosslinkers(Crosslinker.LoadCrosslinkers(crosslinkerLocation));

            // load custom crosslinkers
            string customCrosslinkerLocation = Path.Combine(DataDir, @"Data", @"CustomCrosslinkers.tsv");
            if (File.Exists(customCrosslinkerLocation))
            {
                AddCrosslinkers(Crosslinker.LoadCrosslinkers(customCrosslinkerLocation));
            }

            OGlycanLocations = new List<string>();
            foreach (var glycanFile in Directory.GetFiles(Path.Combine(DataDir, @"Data", @"OGlycan")))
            {
                OGlycanLocations.Add(glycanFile);
            }
            NGlycanLocations = new List<string>();
            foreach (var glycanFile in Directory.GetFiles(Path.Combine(DataDir, @"Data", @"NGlycan")))
            {
                NGlycanLocations.Add(glycanFile);
            }

            ExperimentalDesignFileName = "ExperimentalDesign.tsv";

            UnimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(Path.Combine(DataDir, @"Data", @"unimod.xml")).ToList();
            PsiModDeserialized = UsefulProteomicsDatabases.Loaders.LoadPsiMod(Path.Combine(DataDir, @"Data", @"PSI-MOD.obo.xml"));
            var formalChargesDictionary = UsefulProteomicsDatabases.Loaders.GetFormalChargesDictionary(PsiModDeserialized);
            UniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(Path.Combine(DataDir, @"Data", @"ptmlist.txt"), formalChargesDictionary).ToList();

            foreach (var modFile in Directory.GetFiles(Path.Combine(DataDir, @"Mods")))
            {
                AddMods(UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(modFile, out var errorMods), false);
            }

            AddMods(UniprotDeseralized.OfType<Modification>(), false);
            AddMods(UnimodDeserialized.OfType<Modification>(), false);

            // populate dictionaries of known mods/proteins for deserialization
            AllModsKnownDictionary = new Dictionary<string, Modification>();
            foreach (Modification mod in AllModsKnown)
            {
                if (!AllModsKnownDictionary.ContainsKey(mod.IdWithMotif))
                {
                    AllModsKnownDictionary.Add(mod.IdWithMotif, mod);
                }
                // no error thrown if multiple mods with this ID are present - just pick one
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

            RefreshAminoAcidDictionary();

            string settingsPath = Path.Combine(DataDir, @"settings.toml");
            if (!File.Exists(settingsPath))
            {
                Toml.WriteFile<GlobalSettings>(new GlobalSettings(), settingsPath);
            }

            GlobalSettings = Toml.ReadFile<GlobalSettings>(settingsPath);
            AllSupportedDissociationTypes = new Dictionary<string, DissociationType> {
                { DissociationType.CID.ToString(), DissociationType.CID },
                { DissociationType.ECD.ToString(), DissociationType.ECD },
                { DissociationType.ETD.ToString(), DissociationType.ETD },
                { DissociationType.HCD.ToString(), DissociationType.HCD },
                { DissociationType.EThcD.ToString(), DissociationType.EThcD },
                { DissociationType.Custom.ToString(), DissociationType.Custom },
                { DissociationType.LowCID.ToString(), DissociationType.LowCID}

                // TODO: allow reading from scan header (autodetect dissociation type)
            };
        }

        public static List<string> ErrorsReadingMods = new List<string>();

        // File locations
        public static string DataDir { get; }

        public static bool StopLoops { get; set; }
        public static string ElementsLocation { get; }
        public static string MetaMorpheusVersion { get; }
        public static GlobalSettings GlobalSettings { get; set; }
        public static IEnumerable<Modification> UnimodDeserialized { get; }
        public static IEnumerable<Modification> UniprotDeseralized { get; }
        public static UsefulProteomicsDatabases.Generated.obo PsiModDeserialized { get; }
        public static IEnumerable<Modification> AllModsKnown { get { return _AllModsKnown.AsEnumerable(); } }
        public static IEnumerable<string> AllModTypesKnown { get { return _AllModTypesKnown.AsEnumerable(); } }
        public static Dictionary<string, Modification> AllModsKnownDictionary { get; set; }
        public static Dictionary<string, DissociationType> AllSupportedDissociationTypes { get; private set; }
        public static List<string> SeparationTypes { get { return _SeparationTypes; } }

        public static string ExperimentalDesignFileName { get; }
        public static IEnumerable<Crosslinker> Crosslinkers { get { return _KnownCrosslinkers.AsEnumerable(); } }
        public static IEnumerable<char> InvalidAminoAcids { get { return _InvalidAminoAcids.AsEnumerable(); } }
        public static List<string> OGlycanLocations { get; }
        public static List<string> NGlycanLocations { get; }

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

        public static void AddSeparationTypes(List<string> separationTypes)
        {
            _SeparationTypes.AddRange(separationTypes);
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

        public static void RefreshAminoAcidDictionary()
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
        public static void StartProcess(string path)
        {
            var p = new Process();
            p.StartInfo = new ProcessStartInfo(path)
            {
                UseShellExecute = true
            };
            p.Start();
        }
    }
}