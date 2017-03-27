using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
using System.Reflection;

namespace EngineLayer
{
    public class GlobalEngineLevelSettings
    {

        #region Private Fields

        private static readonly string elementsLocation = Path.Combine(@"Data", @"elements.dat");
        private static readonly string unimodLocation = Path.Combine(@"Data", @"unimod.xml");
        private static readonly string uniprotLocation = Path.Combine(@"Mods", @"ptmlist.txt");

        #endregion Private Fields

        #region Public Constructors

        static GlobalEngineLevelSettings()
        {
            try
            {
                UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);
                UnimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation).ToList();
                UniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation).ToList();
            }
            catch (WebException)
            {
            }

            MetaMorpheusVersion = Assembly.GetExecutingAssembly().GetName().Version.ToString();
        }

        #endregion Public Constructors

        #region Public Properties

        public static string MetaMorpheusVersion { get; private set; }
        public static IEnumerable<Modification> UnimodDeserialized { get; private set; }
        public static IEnumerable<Modification> UniprotDeseralized { get; private set; }

        #endregion Public Properties

    }
}