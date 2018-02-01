using System.Collections.Generic;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class SearchParameters
    {
        #region Public Constructors

        public SearchParameters()
        {
            DisposeOfFileWhenDone = true;
            AddCompIons = false;
            DoParsimony = true;
            NoOneHitWonders = false;
            ModPeptidesAreUnique = true;
            DoQuantification = false;
            QuantifyPpmTol = 5;
            SearchTarget = true;
            DecoyType = DecoyType.Reverse;
            DoHistogramAnalysis = false;
            HistogramBinTolInDaltons = 0.003;
            DoLocalizationAnalysis = true;
            WritePrunedDatabase = false;
            KeepAllUniprotMods = true;
            MassDiffAcceptorType = MassDiffAcceptorType.OneMM;
            MaxFragmentSize = 30000.0;
            ModsToWriteSelection = new Dictionary<string, int>
            {
                {"ProteinTermMod", 3},
                {"UniProt", 2},
            };
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
        public DecoyType DecoyType { get; set; }
        public MassDiffAcceptorType MassDiffAcceptorType { get; set; }
        public bool WritePrunedDatabase { get; set; }
        public bool KeepAllUniprotMods { get; set; }
        public bool DoLocalizationAnalysis { get; set; }
        public bool DoQuantification { get; set; }
        public SearchType SearchType { get; set; }
        public string CustomMdac { get; set; }
        public double MaxFragmentSize { get; set; }
        public double HistogramBinTolInDaltons { get; set; }
        public Dictionary<string, int> ModsToWriteSelection { get; set; }
        public double MaximumMassThatFragmentIonScoreIsDoubled { get; set; }

        #endregion Public Properties
    }
}