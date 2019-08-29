using System.Collections.Generic;
using UsefulProteomicsDatabases;
using EngineLayer;
using Proteomics;

namespace TaskLayer
{
    public class SearchParameters
    {
        public SearchParameters()
        {
            // default search task parameters
            DisposeOfFileWhenDone = true;
            DoParsimony = true;
            NoOneHitWonders = false;
            ModPeptidesAreDifferent = false;
            DoQuantification = true;
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
            WriteMzId = true;
            WritePepXml = false;

            ModsToWriteSelection = new Dictionary<string, int>
            {
                //Key is modification type.

                //Value is integer 0, 1, 2 and 3 interpreted as:
                //   0:   Do not Write
                //   1:   Write if in DB and Observed
                //   2:   Write if in DB
                //   3:   Write if Observed

                {"N-linked glycosylation", 3},
                {"O-linked glycosylation", 3},
                {"Other glycosylation", 3},
                {"Common Biological", 3},
                {"Less Common", 3},
                {"Metal", 3},
                {"2+ nucleotide substitution", 3},
                {"1 nucleotide substitution", 3},
                {"UniProt", 2},
            };

            WriteDecoys = true;
            WriteContaminants = true;
            LocalFdrCategories = new List<FdrCategory> { FdrCategory.FullySpecific };
        }

        public bool DisposeOfFileWhenDone { get; set; }
        public bool DoParsimony { get; set; }
        public bool ModPeptidesAreDifferent { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool MatchBetweenRuns { get; set; }
        public bool Normalize { get; set; }
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
        public List<FdrCategory> LocalFdrCategories { get; set; }
        public string CustomMdac { get; set; }
        public double MaxFragmentSize { get; set; }
        public double HistogramBinTolInDaltons { get; set; }
        public Dictionary<string, int> ModsToWriteSelection { get; set; }
        public double MaximumMassThatFragmentIonScoreIsDoubled { get; set; }
        public bool WriteMzId { get; set; }
        public bool WritePepXml { get; set; }
        public bool WriteDecoys { get; set; }
        public bool WriteContaminants { get; set; }
        public List<SilacLabel> SilacLabels { get; set; }
    }
}