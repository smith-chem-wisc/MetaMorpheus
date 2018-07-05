using System.Collections.Generic;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    /// <summary>
    /// Search parameters with default settings
    /// </summary>
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
            NumDecoyDatabases = 1;
            ModsToWriteSelection = new Dictionary<string, int>
            {
                //Key is modification type.
                //Select from: Common Fixed; Common Variable; Artifact; Biological; Crosslink; Detached; fallOffN; fallOffC; N-linked glycosylation;
                //O -linked glycosylation; Other glycosylation; missing; Deprecated_Mod; Deprecated_PeptideTermMod; Deprecated_Metal;
                //Deprecated_ProteinTermMod; Deprecated_TrypsinDigestedMod; Deprecated_AnpN_DigestedMode; RNA; 1 nucleotide substitution;
                //2 + nucleotide substitution; Surfactant; TandemMassTag; Unimod; UniProt

                //Value is integer 0, 1, 2 and 3 interpreted as:
                //   0:   Do not Write
                //   1:   Write if in DB and Observed
                //   2:   Write if in DB
                //   3:   Write if Observed

                {"Biological", 3},
                //{"ProteinTermMod", 3},
                {"UniProt", 3},
            };
            OutMzId = true;
            OutPepXML = false;
        }

        public bool DisposeOfFileWhenDone { get; set; }
        public bool DoParsimony { get; set; }
        public bool ModPeptidesAreDifferent { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool MatchBetweenRuns { get; set; }
        public bool Normalize { get; set; }
        public double QuantifyPpmTol { get; set; }
        public bool DoHistogramAnalysis { get; set; }
        public bool SearchTarget { get; set; } = true;
        public DecoyType DecoyType { get; set; }
        public int NumDecoyDatabases { get; set; }
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
        public bool OutMzId { get; set; } = true;
        public bool OutPepXML { get; set; }
    }
}