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
            DoBayesianProteinQuant = false;
            BayesianFoldChangeCutoff = 0.1;
            BayesianRandomSeed = 42;
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
            WriteSdrf = false;
            IterativePepTraining = true;

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
        /// else's uptime unless the user asked for it. Turning it on changes q-values: the
        /// spectral angle is a PEP feature, which a search without a spectral library otherwise
        /// trains at the -1 sentinel for every PSM.
        ///
        /// Applies only to classic and modern peptide searches with HCD or CID fragmentation, since the
        /// model is Prosit 2020 HCD. Semi- and non-specific searches (which compute FDR before
        /// post-search analysis), other dissociation types, and oligo or proteoform searches are
        /// skipped with a warning and a line in results.txt.
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

        /// <summary>
        /// Write an SDRF-Proteomics file describing this experiment alongside the results.
        ///
        /// OPT-IN, and deliberately so. An SDRF's sample half -- organism part, disease, cell type,
        /// replicate structure -- is knowledge no search has; only a human does. A run that emitted
        /// one unconditionally would have to write "not available" wherever it could not find a
        /// value, and a corpus of those passes every validator, produces no drift findings, and
        /// cannot be mined. Opting in is the user saying the sample metadata exists, which is what
        /// makes it worth naming every gap before the run starts.
        /// </summary>
        public bool WriteSdrf { get; set; }

        /// <summary>
        /// Retrain the PEP model on its own output until the count of accepted target peptides stops growing
        /// (semi-supervised, as in Percolator and mokapot). On by default, after entrapment checks on two datasets.
        /// Off, PEP trains once, on labels from the search-score q-value.
        /// </summary>
        public bool IterativePepTraining { get; set; }

        /// <summary>
        /// Run FlashLFQ's Bayesian protein fold-change analysis after label-free quantification, comparing every
        /// condition in the experimental design against <see cref="BayesianControlCondition"/>, and write
        /// BayesianFoldChangeAnalysis.tsv. Needs a design with at least two conditions; off for SILAC.
        /// Slow: it runs Markov chain Monte Carlo (1000 burn-in + 3000 steps) per protein for every condition
        /// compared, with two model fits per unpaired comparison, so on a full proteome with several conditions it
        /// can add minutes to hours. Its FDR column is estimated over all of FlashLFQ's protein groups, before
        /// the file drops contaminants, high-q groups and UNDEFINED, and is not re-estimated on the rows written.
        /// Set in the task toml only; there is no GUI control yet.
        /// </summary>
        public bool DoBayesianProteinQuant { get; set; }

        /// <summary>
        /// The experimental-design condition every other condition is compared against. Must match one of the
        /// design's conditions exactly, or the Bayesian step is skipped with a warning.
        /// </summary>
        public string BayesianControlCondition { get; set; }

        /// <summary>
        /// The fold change, as a log2 value, below which a protein counts as unchanged when the Bayesian step
        /// estimates the probability of a real change. FlashLFQ's default, 0.1.
        /// </summary>
        public double BayesianFoldChangeCutoff { get; set; }

        /// <summary>
        /// Seed for the Bayesian step's Markov chain Monte Carlo sampling, so that the same input gives the same
        /// fold changes. Fixed at 42 by default.
        /// </summary>
        public int BayesianRandomSeed { get; set; }

        /// <summary>
        /// The ProteomeXchange accession (PXD######) of the public dataset this search re-analyses,
        /// written into the SDRF as comment[proteomexchange accession number]. Null for a search of
        /// data that has not been deposited.
        ///
        /// It is the join key for pooling: without it, a reanalysis SDRF describes a search but not
        /// which experiment it searched, so it cannot be matched back to the deposition or to other
        /// reanalyses of the same data. Supplied rather than inferred -- nothing in a spectra file
        /// names the dataset it was deposited under.
        /// </summary>
        public string ProteomeXchangeAccession { get; set; }
    }
}