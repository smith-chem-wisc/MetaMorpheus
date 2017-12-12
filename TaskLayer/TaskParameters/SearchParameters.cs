using EngineLayer;
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
            ModTypeList = defaultModTypes;
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
        public Dictionary<string, int> ModTypeList { get; set; }

        #endregion Public Properties

        private Dictionary<string, int> defaultModTypes = new Dictionary<string, int>
        {
           {"Common Fixed", 0},
           {"Common Variable", 0 },
           {"Detached",0 },
           {"Glycan",0 },
           {"missing",0 },
           {"Mod",0 },
           {"ProteinTermMod",0 },
           {"PeptideTermMod",0 },
           {"Metal",0 },
           {"TrypsinDigestedMod",1 },
           {"1 nucleotide substitution",1 },
           {"2+ nucleotide substitution",1 },
           {"Surfactant",1 },
           {"Tandem Mass Tag",1 },
           {"Unimod",1 },
           {"UniProt",1 },

        };

    }
}