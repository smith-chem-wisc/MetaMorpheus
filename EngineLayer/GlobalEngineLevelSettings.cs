using Nett;
using Newtonsoft.Json.Linq;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net.Http;

namespace EngineLayer
{
    public static class GlobalEngineLevelSettings
    {
        #region Public Fields

        public static readonly string elementsLocation;
        public static readonly string modsLocation;
        public static readonly string settingsTomlLocation;

        #endregion Public Fields

        #region Private Fields

        private static readonly string unimodLocation;
        private static readonly string uniprotLocation;
        private static readonly string psiModLocation;
        private static readonly string proteasesLocation;

        #endregion Private Fields

        #region Public Constructors

        static GlobalEngineLevelSettings()
        {
            string dir;
            var pathToProgramFiles = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
            if (!String.IsNullOrWhiteSpace(pathToProgramFiles) && Directory.GetCurrentDirectory().Contains(pathToProgramFiles))
            {
                ByInstaller = true;
                dir = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "MetaMorpheus");
            }
            else
            {
                ByInstaller = false;
                dir = AppDomain.CurrentDomain.BaseDirectory;
            }
            elementsLocation = Path.Combine(dir, @"Data", @"elements.dat");
            unimodLocation = Path.Combine(dir, @"Data", @"unimod.xml");
            uniprotLocation = Path.Combine(dir, @"Data", @"ptmlist.txt");
            psiModLocation = Path.Combine(dir, @"Data", @"PSI-MOD.obo.xml");
            settingsTomlLocation = Path.Combine(dir, @"settings.toml");
            proteasesLocation = Path.Combine(dir, @"Data", "proteases.tsv");
            modsLocation = Path.Combine(dir, @"Mods");

            // No exceptions to be caught here, since we are just reading these files!
            UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);
            UnimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation).ToList();
            PsiModDeserialized = UsefulProteomicsDatabases.Loaders.LoadPsiMod(psiModLocation);
            var formalChargesDictionary = UsefulProteomicsDatabases.Loaders.GetFormalChargesDictionary(PsiModDeserialized);
            UniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation, formalChargesDictionary).ToList();

            MetaMorpheusVersion = typeof(GlobalEngineLevelSettings).Assembly.GetName().Version.ToString();
            if (MetaMorpheusVersion.Equals("1.0.0.0"))
            {
                AskAboutUpdating = false;
#if DEBUG
                MetaMorpheusVersion = "Not a release version. DEBUG.";
#else
                MetaMorpheusVersion = "Not a release version.";
#endif
            }
            else
                AskAboutUpdating = Toml.ReadFile(settingsTomlLocation).Get<bool>("AskAboutUpdating");

            ProteaseDictionary = LoadProteaseDictionary();
            AllModsKnown = new List<Modification>();
        }

        #endregion Public Constructors

        #region Public Properties

        public static bool AskAboutUpdating { get; }

        public static bool ByInstaller { get; }

        public static string MetaMorpheusVersion { get; }

        public static IEnumerable<Modification> UnimodDeserialized { get; }

        public static IEnumerable<Modification> UniprotDeseralized { get; }

        public static UsefulProteomicsDatabases.Generated.obo PsiModDeserialized { get; }

        public static Dictionary<string, Protease> ProteaseDictionary { get; }

        public static List<Modification> AllModsKnown { get; }

        public static string NewestVersion { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public static void GetVersionNumbersFromWeb()
        {
            // Attempt to get current MetaMorpheus version
            using (var client = new HttpClient())
            {
                client.DefaultRequestHeaders.Add("User-Agent", "Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.2; WOW64; Trident/6.0)");

                using (var response = client.GetAsync("https://api.github.com/repos/smith-chem-wisc/MetaMorpheus/releases/latest").Result)
                {
                    var json = response.Content.ReadAsStringAsync().Result;
                    JObject deserialized = JObject.Parse(json);
                    var assets = deserialized["assets"].Select(b => b["name"].ToString()).ToList();
                    if (!assets.Contains("MetaMorpheusInstaller.msi") || !assets.Contains("MetaMorpheusGuiDotNetFrameworkAppveyor.zip"))
                        throw new MetaMorpheusException("Necessary files do not exist!");
                    NewestVersion = deserialized["tag_name"].ToString();
                }
            }
        }

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
                    AllModsKnown.Add(ye);
            }
        }

        #endregion Public Methods

        #region Private Methods

        private static Dictionary<string, Protease> LoadProteaseDictionary()
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