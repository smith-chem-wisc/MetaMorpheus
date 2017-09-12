using System.Linq;

namespace EngineLayer
{
    public class SearchParameters
    {
        #region Public Constructors

        public SearchParameters()
        {
            DisposeOfFileWhenDone = true;
            AddCompIons = false;
            DoParsimony = false;
            NoOneHitWonders = false;
            ModPeptidesAreUnique = true;
            DoQuantification = false;
            QuantifyPpmTol = 5;
            SearchTarget = true;
            SearchDecoy = true;
            DoHistogramAnalysis = false;
            DoLocalizationAnalysis = true;
            WritePrunedDatabase = false;
            KeepAllUniprotMods = true;
            MassDiffAcceptor = GlobalEngineLevelSettings.SearchModesKnown.First();
        }

        #endregion Public Constructors

        #region Public Properties

        public bool DisposeOfFileWhenDone { get; set; }
        public bool AddCompIons { get; set; }
        public bool DoParsimony { get; set; }
        public bool ModPeptidesAreUnique { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool MatchBetweenRuns { get; set; }
        public double QuantifyPpmTol { get; set; }
        public bool DoHistogramAnalysis { get; set; }
        public bool SearchTarget { get; set; }
        public bool SearchDecoy { get; set; }
        public MassDiffAcceptor MassDiffAcceptor { get; set; }
        public bool WritePrunedDatabase { get; set; }
        public bool KeepAllUniprotMods { get; set; }
        public bool DoLocalizationAnalysis { get; set; }
        public bool DoQuantification { get; set; }
        public SearchType SearchType { get; set; }

        #endregion Public Properties
    }
}