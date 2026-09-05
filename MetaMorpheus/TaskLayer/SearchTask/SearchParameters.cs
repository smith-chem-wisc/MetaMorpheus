using System.Collections.Generic;
using UsefulProteomicsDatabases;
using EngineLayer;
using Omics.Modifications;

namespace TaskLayer
{
    public class SearchParameters
    {
        /// <summary>
        /// Default maximum fragment size in Daltons used for indexing. This value is shared across
        /// all task types that require fragment indexing (Search, Calibration, CrossLink, Glyco).
        /// </summary>
        public const double DefaultMaxFragmentSize = 30000.0;

        public SearchParameters()
        {
            // default search task parameters
            DisposeOfFileWhenDone = true;
            DoParsimony = true;
            NoOneHitWonders = false;
            ModPeptidesAreDifferent = false;
            DoLabelFreeQuantification = true;
            UseSharedPeptidesForLFQ = false;
            QuantifyPpmTol = 5;
            MbrFdrThreshold = 0.01;
            SearchTarget = true;
            DecoyType = DecoyType.Reverse;
            DoHistogramAnalysis = false;
            HistogramBinTolInDaltons = 0.003;
            DoLocalizationAnalysis = true;
            WritePrunedDatabase = false;
            KeepAllUniprotMods = true;
            MassDiffAcceptorType = MassDiffAcceptorType.OneMM;
            MaxFragmentSize = DefaultMaxFragmentSize;
            MinAllowedInternalFragmentLength = 0;
            UsePredictedSpectraForSpectralAngle = false;
            WriteMzId = true;
            WritePepXml = false;
            IncludeModMotifInMzid = false;
            WriteDigestionProductCountFile = false;
            WriteTargetDecoyFasta = false;

            ModsToWriteSelection = DefaultModsToWriteSelection();

            WriteHighQValuePsms = true;
            WriteDecoys = true;
            WriteContaminants = true;
            WriteIndividualFiles = true;
            LocalFdrCategories = new List<FdrCategory> { FdrCategory.FullySpecific };
            TCAmbiguity = TargetContaminantAmbiguity.RemoveContaminant;
        }

        public bool DisposeOfFileWhenDone { get; set; }
        public bool DoParsimony { get; set; }
        public bool ModPeptidesAreDifferent { get; set; }
        public bool NoOneHitWonders { get; set; }
        public bool MatchBetweenRuns { get; set; }
        public double MbrFdrThreshold { get; set; }
        public bool Normalize { get; set; }
        public double QuantifyPpmTol { get; set; }
        public bool DoHistogramAnalysis { get; set; }
        public bool SearchTarget { get; set; }
        public DecoyType DecoyType { get; set; }
        public MassDiffAcceptorType MassDiffAcceptorType { get; set; }
        public bool WritePrunedDatabase { get; set; }
        public bool KeepAllUniprotMods { get; set; }
        public bool DoLocalizationAnalysis { get; set; }
        public bool DoLabelFreeQuantification { get; set; }
        public bool UseSharedPeptidesForLFQ { get; set; }
        public bool DoMultiplexQuantification { get; set; }
        public string MultiplexModId { get; set; }
        public SearchType SearchType { get; set; }
        public List<FdrCategory> LocalFdrCategories { get; set; }
        public string CustomMdac { get; set; }
        public double MaxFragmentSize { get; set; }
        public int MinAllowedInternalFragmentLength { get; set; } //0 means "no internal fragments"
        public double HistogramBinTolInDaltons { get; set; }

        /// <summary>
        /// The default modification types written to a pruned database, keyed by modification type.
        /// Values are 0 do not write, 1 write if in the database and observed, 2 write if in the database,
        /// 3 write if observed. A fresh dictionary each call, since callers mutate their own copy.
        /// </summary>
        /// <remarks>
        /// Shared with <see cref="GlycoSearchParameters"/>, which needs the same protein defaults.
        /// <see cref="RnaSearchParameters"/> deliberately replaces it with an RNA-specific set.
        /// </remarks>
        public static Dictionary<string, int> DefaultModsToWriteSelection() => new Dictionary<string, int>
        {
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
        public Dictionary<string, int> ModsToWriteSelection { get; set; }
        public double MaximumMassThatFragmentIonScoreIsDoubled { get; set; }
        public bool WriteMzId { get; set; }
        public bool WritePepXml { get; set; }
        public bool WriteHighQValuePsms { get; set; }
        public bool WriteDecoys { get; set; }
        public bool WriteContaminants { get; set; }
        public bool WriteIndividualFiles { get; set; }
        public bool WriteSpectralLibrary { get; set; }
        /// <summary>
        /// Opt in to filling missing spectral angles with Prosit-predicted spectra. Off by
        /// default because it is a call to a third-party web service (Koina) on every search:
        /// a search that would otherwise run offline should not start depending on someone
        /// else's uptime unless the user asked for it. Angles from a real spectral library are
        /// unaffected and still computed either way.
        /// </summary>
        public bool UsePredictedSpectraForSpectralAngle { get; set; }
        public bool UpdateSpectralLibrary { get; set; }
        public bool CompressIndividualFiles { get; set; }
        public List<SilacLabel> SilacLabels { get; set; }
        public SilacLabel StartTurnoverLabel { get; set; } //used for SILAC turnover experiments
        public SilacLabel EndTurnoverLabel { get; set; } //used for SILAC turnover experiments
        public TargetContaminantAmbiguity TCAmbiguity { get; set; }
        public bool IncludeModMotifInMzid { get; set; }
        public bool WriteDigestionProductCountFile { get; set; }
        public bool WriteTargetDecoyFasta { get; set; }
    }
}