using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer
{
    public static class GlobalVariables
    {
        #region Private Fields

        private static List<Modification> allModsKnown;

        #endregion Private Fields

        #region Public Constructors

        static GlobalVariables()
        {
            #region Determine MetaMorpheusVersion

            MetaMorpheusVersion = typeof(GlobalVariables).Assembly.GetName().Version.ToString();

            if (MetaMorpheusVersion.Equals("1.0.0.0"))
            {
#if DEBUG
                MetaMorpheusVersion = "Not a release version. DEBUG.";
#else
                MetaMorpheusVersion = "Not a release version.";
#endif
            }

            #endregion Determine MetaMorpheusVersion

            #region Figure out DataDir

            {
                var pathToProgramFiles = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
                if (!String.IsNullOrWhiteSpace(pathToProgramFiles) && AppDomain.CurrentDomain.BaseDirectory.Contains(pathToProgramFiles) && !AppDomain.CurrentDomain.BaseDirectory.Contains("Jenkins"))
                    DataDir = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "MetaMorpheus");
                else
                    DataDir = AppDomain.CurrentDomain.BaseDirectory;
            }

            #endregion Figure out DataDir

            ElementsLocation = Path.Combine(DataDir, @"Data", @"elements.dat");
            UsefulProteomicsDatabases.Loaders.LoadElements(ElementsLocation);

            UnimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(Path.Combine(DataDir, @"Data", @"unimod.xml")).ToList();
            PsiModDeserialized = UsefulProteomicsDatabases.Loaders.LoadPsiMod(Path.Combine(DataDir, @"Data", @"PSI-MOD.obo.xml"));
            var formalChargesDictionary = UsefulProteomicsDatabases.Loaders.GetFormalChargesDictionary(PsiModDeserialized);
            UniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(Path.Combine(DataDir, @"Data", @"ptmlist.txt"), formalChargesDictionary).ToList();

            allModsKnown = new List<Modification>();
            foreach (var modFile in Directory.GetFiles(Path.Combine(DataDir, @"Mods")))
                AddMods(UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(modFile));
            AddMods(UnimodDeserialized.OfType<ModificationWithLocation>());
            AddMods(UniprotDeseralized.OfType<ModificationWithLocation>());

            GlobalSettings = Toml.ReadFile<GlobalSettings>(Path.Combine(DataDir, @"settings.toml"));

            ProteaseDictionary = LoadProteaseDictionary(Path.Combine(DataDir, @"Data", "proteases.tsv"));
        }

        #endregion Public Constructors

        #region Public Properties

        // File locations
        public static string DataDir { get; }

        public static string ElementsLocation { get; }
        public static string MetaMorpheusVersion { get; }
        public static IGlobalSettings GlobalSettings { get; }
        public static IEnumerable<Modification> UnimodDeserialized { get; }
        public static IEnumerable<Modification> UniprotDeseralized { get; }
        public static UsefulProteomicsDatabases.Generated.obo PsiModDeserialized { get; }
        public static IReadOnlyDictionary<string, Protease> ProteaseDictionary { get; }
        public static IEnumerable<Modification> AllModsKnown { get { return allModsKnown.AsEnumerable(); } }

        #endregion Public Properties

        #region Public Methods

        public static void AddMods(IEnumerable<Modification> enumerable)
        {
            foreach (var ye in enumerable)
            {
                if (string.IsNullOrEmpty(ye.modificationType) || string.IsNullOrEmpty(ye.id))
                    throw new MetaMorpheusException(ye.ToString() + Environment.NewLine + " has null or empty modification type");
                if (AllModsKnown.Any(b => b.id.Equals(ye.id) && b.modificationType.Equals(ye.modificationType) && !b.Equals(ye)))
                    throw new MetaMorpheusException("Modification id and type are equal, but some fields are not! Please modify/remove one of the modifications: " + Environment.NewLine + Environment.NewLine + ye.ToString() + Environment.NewLine + Environment.NewLine + " has same and id and modification type as " + Environment.NewLine + Environment.NewLine + AllModsKnown.First(b => b.id.Equals(ye.id) && b.modificationType.Equals(ye.modificationType)) + Environment.NewLine + Environment.NewLine);
                else if (AllModsKnown.Any(b => b.id.Equals(ye.id) && b.modificationType.Equals(ye.modificationType)))
                    continue;
                else
                    allModsKnown.Add(ye);
            }
        }

        #endregion Public Methods

        #region Private Methods

        private static Dictionary<string, Protease> LoadProteaseDictionary(string proteasesLocation)
        {
            Dictionary<string, Protease> dict = new Dictionary<string, Protease>();
            using (StreamReader proteases = new StreamReader(proteasesLocation))
            {
                proteases.ReadLine();

                while (proteases.Peek() != -1)
                {
                    string line = proteases.ReadLine();
                    string[] fields = line.Split('\t');

                    string name = fields[0];
                    string[] sequences_inducing_cleavage = fields[1].Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                    string[] sequences_preventing_cleavage = fields[2].Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                    var cleavage_terminus = (TerminusType)Enum.Parse(typeof(TerminusType), fields[3], true);
                    var cleavage_specificity = (CleavageSpecificity)Enum.Parse(typeof(CleavageSpecificity), fields[4], true);
                    string psi_ms_accession_number = fields[5];
                    string psi_ms_name = fields[6];
                    string site_regexp = fields[7];
                    var protease = new Protease(name, sequences_inducing_cleavage, sequences_preventing_cleavage, cleavage_terminus, cleavage_specificity, psi_ms_accession_number, psi_ms_name, site_regexp);
                    dict.Add(protease.Name, protease);
                }
            }
            return dict;
        }

        #endregion Private Methods
    }
}