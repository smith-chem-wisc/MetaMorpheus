using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
using System.Reflection;

namespace EngineLayer
{
    public class GlobalEngineLevelSettings
    {
        #region Public Fields

        public static readonly string elementsLocation = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"Data", @"elements.dat");

        #endregion Public Fields

        #region Private Fields

        private static readonly string unimodLocation = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"Data", @"unimod.xml");
        private static readonly string uniprotLocation = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"Data", @"ptmlist.txt");
        private static readonly string psiModLocation = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"Data", @"PSI-MOD.obo.xml");

        #endregion Private Fields

        #region Public Constructors

        static GlobalEngineLevelSettings()
        {
            try
            {
                UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);
                UnimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation).ToList();
                PsiModDeserialized = UsefulProteomicsDatabases.Loaders.LoadPsiMod(psiModLocation);
                var formalChargesDictionary = UsefulProteomicsDatabases.Loaders.GetFormalChargesDictionary(PsiModDeserialized);
                UniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation, formalChargesDictionary).ToList();
            }
            catch (WebException)
            {
            }

            MetaMorpheusVersion = Assembly.GetExecutingAssembly().GetName().Version.ToString();
            if (MetaMorpheusVersion.Equals("1.0.0.0"))
                MetaMorpheusVersion = "Not a release version";
        }

        #endregion Public Constructors

        #region Public Properties

        public static string MetaMorpheusVersion { get; }
        public static IEnumerable<Modification> UnimodDeserialized { get; }
        public static IEnumerable<Modification> UniprotDeseralized { get; }
        public static UsefulProteomicsDatabases.Generated.obo PsiModDeserialized { get; }

        #endregion Public Properties
    }
}