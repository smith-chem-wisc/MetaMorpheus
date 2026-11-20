using System.Collections.Generic;
using UsefulProteomicsDatabases;
using EngineLayer;
using EngineLayer.Localization;
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
            WriteSdrf = false;
            IterativePepTraining = true;

            ModsToWriteSelection = DefaultModsToWriteSelection();

            WriteHighQValuePsms = true;
            WriteDecoys = true;
            WriteContaminants = true;
            WriteIndividualFiles = false;
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

        /// <summary>
        /// Estimate, for each match, the probability of each placement of its modifications (and, in an open search, of
        /// its mass difference), and a localization q-value: the false localization rate of the confident matches as
        /// sure of their placement as it is (PTMiner, An et al., MCP 2019). Off by default: it reads every spectra file
        /// a second time and keeps each match's matched peaks until the end of the search.
        /// </summary>
        public bool DoLocalizationProbabilities { get; set; }

        /// <summary>
        /// With <see cref="DoLocalizationProbabilities"/>, which residues an identified modification may be placed on:
        /// those its motif allows (the default), or any residue, for experiments where its chemistry is unknown.
        /// </summary>
        public LocalizationSiteRule LocalizationSiteRule { get; set; }

        /// <summary>
        /// With <see cref="DoLocalizationProbabilities"/>, match ion-trap scans at the low-resolution product tolerance
        /// (<see cref="EngineLayer.CommonParameters.ProductMassTolerance_LowRes"/>) when collecting localization evidence,
        /// and learn their fragment statistics apart from other scan types. For methods that put ion-trap and Orbitrap
        /// MS2 in one file. Off by default.
        /// </summary>
        public bool LocalizeIonTrapScansSeparately { get; set; }

        /// <summary>
        /// With <see cref="DoLocalizationProbabilities"/>, learn how often each modification sits on each residue from the
        /// matches (PTMiner's prior) instead of weighing every allowed residue alike. Off by default: a learned prior
        /// helps accuracy but leaves the residues it disfavours overconfident.
        /// </summary>
        public bool LearnLocalizationPrior { get; set; }

        /// <summary>
        /// With <see cref="DoLocalizationProbabilities"/>, the power each placement's likelihood is raised to, tempering
        /// the fragment evidence, which the model treats as independent when it is not. 1 leaves it untempered.
        /// </summary>
        public double LocalizationLikelihoodExponent { get; set; } = LocalizationProbabilityEngine.RecommendedLikelihoodExponent;

        /// <summary>
        /// With <see cref="DoLocalizationProbabilities"/>, the retention-time window, in minutes, within which another
        /// identification of the same peptide with its modifications placed differently marks a match as having a
        /// co-eluting isomer, which then gets no localization q-value. 0 checks nothing.
        /// </summary>
        public double LocalizationCoElutionWindow { get; set; } = LocalizationProbabilityEngine.RecommendedCoElutionWindow;
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
        /// <summary>
        /// On a search of more than one spectra file, write a PSM, peptide, protein-group and quantified-peak
        /// table for each file to "Individual File Results", and each file's counts to the results summary
        /// (results.txt, allResults.txt). The per-file PSM and peptide tables re-run FDR on that file's
        /// matches alone. Off by default; a task TOML that sets it to true still writes them. Per-file .mzID
        /// and .pep.XML are governed by <see cref="WriteMzId"/> and <see cref="WritePepXml"/>, not by this
        /// setting.
        /// </summary>
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

        /// <summary>
        /// An SDRF describing the samples this search's files came from, or null. With <see cref="WriteSdrf"/>, the
        /// SDRF the search writes then copies each searched file's sample columns from it cell for cell (source name,
        /// characteristics, material type, factor values and sample comments), instead of building them from
        /// ExperimentalDesign.tsv or TmtDesign.txt. Statistics reads per-sample facts (individual, age, batch) only
        /// from the written SDRF, and the curator's cells are the only record of them (sdrf D50).
        ///
        /// It describes the SDRF's samples only; it is not the experimental design. Quantification, and the SDRF's
        /// technical replicate, fraction and channels, still come from ExperimentalDesign.tsv or TmtDesign.txt. To make
        /// that design from an SDRF, the command line's --sdrfDesign writes it beside the spectra files; give it the same
        /// SDRF, and the design and the written SDRF's sample columns agree.
        ///
        /// Give a full path: a relative one is read against the working directory, not the task file's folder. Set in
        /// the task file by a pipeline that re-analyses a deposit; the GUI has no field for it.
        /// </summary>
        public string InputSdrfPath { get; set; }
    }
}