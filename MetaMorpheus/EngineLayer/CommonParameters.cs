using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Nett;
using Omics.Digestion;
using Transcriptomics.Digestion;
using EngineLayer.DIA;
using Transcriptomics;
using EngineLayer.FdrAnalysis;
using EngineLayer.SpectrumMatch;

namespace EngineLayer
{
    public class CommonParameters
    {
        // this parameterless constructor needs to exist to read the toml.
        // if you can figure out a way to get rid of it, feel free...
        public CommonParameters()
            : this(digestionParams: null)
        {
            // A task TOML is read into this instance, setting only the keys the file holds. Until it names the
            // identification filter or its threshold, each is "not named" and is resolved from what the file does hold
            // (IdentificationFilter): a file saved before FilterThreshold existed is read by its old thresholds.
            _filterThresholdNamed = false;
        }

        public CommonParameters(
            string taskDescriptor = null, 
            DissociationType dissociationType = DissociationType.HCD, 
            DissociationType ms2childScanDissociationType = DissociationType.Unknown, 
            DissociationType ms3childScanDissociationType = DissociationType.Unknown, 
            string separationType = "HPLC", 
            bool doPrecursorDeconvolution = true,
            bool useProvidedPrecursorInfo = true, 
            double deconvolutionIntensityRatio = 3, 
            int deconvolutionMaxAssumedChargeState = 12, 
            bool reportAllAmbiguity = true,
            bool addCompIons = false, 
            int totalPartitions = 1,
            FilterType? identificationFilter = null,
            double filterThreshold = DefaultFilterThreshold,
            double qValueCutoffForPepCalculation = 0.005,
            double scoreCutoff = 5, 
            int? numberOfPeaksToKeepPerWindow = 200, 
            double? minimumAllowedIntensityRatioToBasePeak = 0.01, 
            double? windowWidthThomsons = null, 
            int? numberOfWindows = null, 
            bool normalizePeaksAccrossAllWindows = false, 
            bool trimMs1Peaks = false,
            bool trimMsMsPeaks = true, 
            Tolerance productMassTolerance = null, 
            Tolerance precursorMassTolerance = null,
            Tolerance productMassTolerance_LowRes = null,
            Tolerance deconvolutionMassTolerance = null,
            int maxThreadsToUsePerFile = -1, 
            IDigestionParams digestionParams = null, 
            IEnumerable<(string, string)> listOfModsVariable = null, 
            IEnumerable<(string, string)> listOfModsFixed = null,
            bool assumeOrphanPeaksAreZ1Fragments = true, 
            int maxHeterozygousVariants = 4, 
            int minVariantDepth = 1, 
            bool addTruncations = false,
            DeconvolutionParameters precursorDeconParams = null,
            DeconvolutionParameters productDeconParams = null,
            bool useMostAbundantPrecursorIntensity = true,
            DIAparameters diaParameters = null,
            IFragmentationParams fragmentationParams = null,
            PrecursorMassMatchMode precursorMassMatchMode = PrecursorMassMatchMode.Monoisotopic,
                string rtPredictorName = RTPredictorNames.Chronologer,
            DoubleRange retentionTimeRange = null,
            string fivePrimeTerminus = null,
            string threePrimeTerminus = null)

        {
            retentionTimeRange ??= new DoubleRange(0, double.MaxValue);
            if (double.IsNaN(retentionTimeRange.Minimum) || retentionTimeRange.Minimum < 0)
                throw new ArgumentOutOfRangeException(nameof(retentionTimeRange), "Minimum retention time must be non-negative.");

            if (double.IsNaN(retentionTimeRange.Maximum) || retentionTimeRange.Maximum < retentionTimeRange.Minimum)
                throw new ArgumentOutOfRangeException(nameof(retentionTimeRange), "Maximum retention time must be greater than or equal to minimum retention time.");

            TaskDescriptor = taskDescriptor;
            DoPrecursorDeconvolution = doPrecursorDeconvolution;
            UseProvidedPrecursorInfo = useProvidedPrecursorInfo;
            DeconvolutionIntensityRatio = deconvolutionIntensityRatio;
            ReportAllAmbiguity = reportAllAmbiguity;
            AddCompIons = addCompIons;
            TotalPartitions = totalPartitions;
            // Null leaves the filter unnamed: the default for the analyte the digestion parameters search (IdentificationFilter)
            if (identificationFilter is FilterType filter)
            {
                IdentificationFilter = filter;
            }
            FilterThreshold = filterThreshold;
            QValueCutoffForPepCalculation = qValueCutoffForPepCalculation;
            ScoreCutoff = scoreCutoff;
            NumberOfPeaksToKeepPerWindow = numberOfPeaksToKeepPerWindow;
            MinimumAllowedIntensityRatioToBasePeak = minimumAllowedIntensityRatioToBasePeak;
            WindowWidthThomsons = windowWidthThomsons;
            NumberOfWindows = numberOfWindows;
            NormalizePeaksAccrossAllWindows = normalizePeaksAccrossAllWindows;
            TrimMs1Peaks = trimMs1Peaks;
            TrimMsMsPeaks = trimMsMsPeaks;
            MaxThreadsToUsePerFile = maxThreadsToUsePerFile == -1 ? Environment.ProcessorCount > 1 ? Environment.ProcessorCount - 1 : 1 : maxThreadsToUsePerFile;
            ProductMassTolerance = productMassTolerance ?? new PpmTolerance(20);
            PrecursorMassTolerance = precursorMassTolerance ?? new PpmTolerance(5);
            ProductMassTolerance_LowRes = productMassTolerance_LowRes ?? new AbsoluteTolerance(0.35);
            DeconvolutionMassTolerance = deconvolutionMassTolerance ?? new PpmTolerance(4);
            DigestionParams = digestionParams ?? new DigestionParams();
            DissociationType = dissociationType;
            SeparationType = separationType;
            MS2ChildScanDissociationType = ms2childScanDissociationType;
            MS3ChildScanDissociationType = ms3childScanDissociationType;
            UseMostAbundantPrecursorIntensity = useMostAbundantPrecursorIntensity;
            PrecursorMassMatchMode = precursorMassMatchMode;
            AssumeOrphanPeaksAreZ1Fragments = assumeOrphanPeaksAreZ1Fragments;
            MaxHeterozygousVariants = maxHeterozygousVariants;
            MinVariantDepth = minVariantDepth;
            AddTruncations = addTruncations;
            DIAparameters = diaParameters;
            RetentionTimeRange = retentionTimeRange;

            // product maximum charge state of 10 is a preexisting hard-coded value in MetaMorpheus
            if (deconvolutionMaxAssumedChargeState > 0) // positive mode
            {
                PrecursorDeconvolutionParameters = precursorDeconParams ?? new ClassicDeconvolutionParameters(1,
                    deconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance.Value, deconvolutionIntensityRatio);
                ProductDeconvolutionParameters = productDeconParams ?? new ClassicDeconvolutionParameters(1,
                    10, DeconvolutionMassTolerance.Value, deconvolutionIntensityRatio);
            }
            else // negative mode
            {
                PrecursorDeconvolutionParameters = precursorDeconParams ?? new ClassicDeconvolutionParameters(deconvolutionMaxAssumedChargeState,
                    -1, DeconvolutionMassTolerance.Value, deconvolutionIntensityRatio, Polarity.Negative);
                ProductDeconvolutionParameters = productDeconParams ?? new ClassicDeconvolutionParameters(-10,
                    -1, DeconvolutionMassTolerance.Value, deconvolutionIntensityRatio, Polarity.Negative);
            }

            if (digestionParams is RnaDigestionParams)
            {
                ListOfModsVariable = listOfModsVariable ?? new List<(string, string)> { ("Digestion Termini", "Cyclic Phosphate on X") };
                ListOfModsFixed = listOfModsFixed ?? new List<(string, string)>();
                PrecursorDeconvolutionParameters.AverageResidueModel = new OxyriboAveragine();
                ProductDeconvolutionParameters.AverageResidueModel = new OxyriboAveragine();
                FragmentationParameters = fragmentationParams ?? RnaFragmentationParams.Default;
            }
            else
            {
                ListOfModsVariable = listOfModsVariable ?? new List<(string, string)> { ("Common Variable", "Oxidation on M") };
                ListOfModsFixed = listOfModsFixed ?? new List<(string, string)> { ("Common Fixed", "Carbamidomethyl on C"), ("Common Fixed", "Carbamidomethyl on U") };
                FragmentationParameters = fragmentationParams ?? new FragmentationParams();
            }

            RTPredictorName = rtPredictorName;

            // Terminal caps only mean something for RNA; null means "use mzLib's default cap".
            bool isRna = digestionParams is RnaDigestionParams;
            FivePrimeTerminus = isRna ? TryParseFormula(fivePrimeTerminus) : null;
            ThreePrimeTerminus = isRna ? TryParseFormula(threePrimeTerminus) : null;

            CustomIons = digestionParams.ProductsFromDissociationType()[DissociationType.Custom];

            // reset custom fragmentation product types to default empty list
            digestionParams.ProductsFromDissociationType()[DissociationType.Custom] = new List<ProductType>() { };
        }

        // Notes:
        // 1) Any new property must not be nullable (such as int?) or else if it is null,
        //    the null setting will not be written to a toml
        //    and the default will override (so it's okay ONLY if the default is null)
        // 2) All setters should be private unless necessary
        // 3) If you add a new property here, you must add it to MetaMorpheusTask.cs/SetAllFileSpecificCommonParams !!
        //    If you forget this, and the user uses file-specific parameters, your settings will be overwritten by default values!

        public string TaskDescriptor { get; private set; }
        public int MaxThreadsToUsePerFile { get; set; }
        public IEnumerable<(string, string)> ListOfModsFixed { get; private set; }
        public IEnumerable<(string, string)> ListOfModsVariable { get; private set; }
        public bool DoPrecursorDeconvolution { get; set; }
        public bool UseProvidedPrecursorInfo { get; set; }
        [TomlIgnore] public double DeconvolutionIntensityRatio { get; private set; }
        public int DeconvolutionMaxAssumedChargeState
        {
            get => PrecursorDeconvolutionParameters.MaxAssumedChargeState;
            private set => PrecursorDeconvolutionParameters.MaxAssumedChargeState = value;
        }
        public DeconvolutionParameters PrecursorDeconvolutionParameters { get; private set; }
        public DeconvolutionParameters ProductDeconvolutionParameters { get; private set; }
        [TomlIgnore] public Tolerance DeconvolutionMassTolerance { get; private set; }
        public int TotalPartitions { get; set; }
        public Tolerance ProductMassTolerance { get; set; } // public setter required for calibration task
        public Tolerance ProductMassTolerance_LowRes { get; set; }// Use a wider mass tolerance for lower-resolution analyzers (e.g., ion traps). For now, this is an independent parameter used only in glyco tasks and is not modified by the calibration task.
        public Tolerance PrecursorMassTolerance { get; set; } // public setter required for calibration task
        public bool AddCompIons { get; set; }

        private FilterType _identificationFilter;
        private bool _identificationFilterNamed;
        private double _filterThreshold;
        private bool _filterThresholdNamed;

        // What a task TOML saved before FilterThreshold existed held instead: recorded by the task TOML reader
        // (RecordLegacyThreshold), never written back. _oldTaskFile: it held at least one of the old keys.
        private double? _legacyQValueThreshold;
        private double? _legacyPepQValueThreshold;
        private bool _oldTaskFile;

        /// <summary>The key a task TOML saved before <see cref="FilterThreshold"/> existed held the q-value threshold under.</summary>
        public const string LegacyQValueThresholdKey = "QValueThreshold";

        /// <summary>The key a task TOML saved before <see cref="FilterThreshold"/> existed held the PEP q-value threshold under.</summary>
        public const string LegacyPepQValueThresholdKey = "PepQValueThreshold";

        /// <summary>
        /// The keys a task TOML saved before QValueThreshold and PepQValueThreshold existed (MetaMorpheus before June 2023,
        /// #2285) held its output filters under. Their values filtered only what was written, and have been ignored since;
        /// they still are. They mark the file as an old task file, which keeps the q-value notch at 0.01 it has searched
        /// with since (<see cref="IdentificationFilter"/>).
        /// </summary>
        public static readonly IReadOnlyList<string> LegacyOutputFilterKeys = new[] { "QValueOutputFilter", "PepQValueOutputFilter" };

        /// <summary>What a task TOML that held only one of the old thresholds meant by QValueThreshold: its old default.</summary>
        public const double LegacyQValueThresholdDefault = 0.01;

        /// <summary>What a task TOML that held only one of the old thresholds meant by PepQValueThreshold: its old default.</summary>
        public const double LegacyPepQValueThresholdDefault = 1.0;

        /// <summary>
        /// The value that decides whether a match or peptide counts as identified, at every step of a search: the
        /// tables, parsimony, protein scoring and FDR, quantification, the spectral library, the pruned database, the
        /// Spectrum Match Count column and every count reported. It is held to <see cref="FilterThreshold"/>. The user
        /// picks one of three:
        /// <list type="bullet">
        /// <item><description><see cref="FilterType.PepQValue"/>: the PEP q-value. A match whose PEP was not trained is
        /// judged by its q-value notch at the same threshold, and the search says so.</description></item>
        /// <item><description><see cref="FilterType.QValueNotch"/>: the q-value counted within each mass-difference
        /// notch.</description></item>
        /// <item><description><see cref="FilterType.QValue"/>: the pooled q-value, counted over every notch together, and
        /// nothing else.</description></item>
        /// </list>
        /// Unnamed, the filter is the default for the analyte the <see cref="DigestionParams"/> search, read when it is asked
        /// for (a task TOML sets its digestion parameters after constructing this instance):
        /// <see cref="DefaultIdentificationFilter(bool, bool)"/>, the PEP q-value bottom-up and the q-value notch for
        /// top-down and RNA.
        /// <para>A task TOML saved before <see cref="FilterThreshold"/> existed holds QValueThreshold and PepQValueThreshold
        /// instead (every file MetaMorpheus wrote from June 2023 holds both), and keeps what it meant: the filter it names
        /// at that filter's old threshold (the PEP q-value at PepQValueThreshold, otherwise QValueThreshold); naming none,
        /// the PEP q-value when PepQValueThreshold is below QValueThreshold (the old GUI's PEP box saved PepQValueThreshold
        /// 0.01 and QValueThreshold 1.0; <see cref="IdentificationFilterFromLegacyThresholds"/>), the q-value notch at
        /// QValueThreshold otherwise. A key the file lacks means its old default (<see cref="LegacyQValueThresholdDefault"/>,
        /// <see cref="LegacyPepQValueThresholdDefault"/>); a file older still holds only
        /// <see cref="LegacyOutputFilterKeys"/>, and so means the q-value notch at 0.01. Saving it writes the filter and
        /// its threshold, so the rule applies once per old file. A filter the TOML names, the constructor is given or the
        /// setter sets, always wins; so does a FilterThreshold the TOML holds.</para>
        /// </summary>
        public FilterType IdentificationFilter
        {
            get => ResolveIdentificationFilter().Filter;
            private set
            {
                _identificationFilter = value;
                _identificationFilterNamed = true;
            }
        }

        /// <summary>
        /// The one threshold <see cref="IdentificationFilter"/> is held to: a match or peptide passes when the filter's value
        /// is at or below it. 0.01 by default. If SearchParameters.WriteHighQValuePsms is false, only the matches and
        /// peptides that pass are written to the .psmtsv files. A task TOML saved before this setting existed is read by
        /// its old thresholds (<see cref="IdentificationFilter"/>).
        /// </summary>
        public double FilterThreshold
        {
            get => ResolveIdentificationFilter().Threshold;
            private set
            {
                _filterThreshold = value;
                _filterThresholdNamed = true;
            }
        }

        /// <summary>The <see cref="FilterThreshold"/> a search gets when it names none.</summary>
        public const double DefaultFilterThreshold = 0.01;

        /// <summary>
        /// The identification filter a search gets when it names none: the q-value notch for RNA
        /// (<paramref name="isRna"/>) and for top-down (<paramref name="isTopDown"/>, the "top-down" protease), the PEP
        /// q-value bottom-up. The task windows start a top-down task from it, as from
        /// <see cref="PepAnalysisEngine.DefaultQValueCutoffForPepCalculation"/>.
        /// </summary>
        public static FilterType DefaultIdentificationFilter(bool isRna, bool isTopDown) =>
            isRna || isTopDown ? FilterType.QValueNotch : FilterType.PepQValue;

        private (FilterType Filter, double Threshold) ResolveIdentificationFilter()
        {
            // FilterThreshold named, or not an old task file: the filter as named, or the analyte's default
            if (_filterThresholdNamed || !_oldTaskFile)
            {
                if (_identificationFilterNamed)
                {
                    return (_identificationFilter, _filterThreshold);
                }
                AnalyteType analyte = DetermineAnalyteType();
                return (DefaultIdentificationFilter(analyte == AnalyteType.Oligo, analyte == AnalyteType.Proteoform), _filterThreshold);
            }

            // A task TOML saved before FilterThreshold existed: what its old thresholds meant
            double qValueThreshold = _legacyQValueThreshold ?? LegacyQValueThresholdDefault;
            double pepQValueThreshold = _legacyPepQValueThreshold ?? LegacyPepQValueThresholdDefault;
            if (_identificationFilterNamed)
            {
                return (_identificationFilter, _identificationFilter == FilterType.PepQValue ? pepQValueThreshold : qValueThreshold);
            }
            return pepQValueThreshold < qValueThreshold
                ? (FilterType.PepQValue, pepQValueThreshold)
                : (FilterType.QValueNotch, qValueThreshold);
        }

        /// <summary>
        /// Records a key a task TOML saved before <see cref="FilterThreshold"/> existed holds, for
        /// <see cref="IdentificationFilter"/> and <see cref="FilterThreshold"/> to read: one of the two old thresholds
        /// (<see cref="LegacyQValueThresholdKey"/>, <see cref="LegacyPepQValueThresholdKey"/>), with its value, or one of the
        /// older output filters (<see cref="LegacyOutputFilterKeys"/>), whose value is ignored. The task TOML reader
        /// (MetaMorpheusTask.tomlConfig) calls it for each key CommonParameters has no property for. Returns false, recording
        /// nothing, for any other key.
        /// </summary>
        /// <exception cref="MetaMorpheusException">The key is one of the two thresholds, and its value is not a number.</exception>
        public bool RecordLegacyThreshold(string key, TomlObject value)
        {
            if (LegacyOutputFilterKeys.Contains(key))
            {
                _oldTaskFile = true;
                return true;
            }
            if (key is not (LegacyQValueThresholdKey or LegacyPepQValueThresholdKey))
            {
                return false;
            }

            double threshold = value switch
            {
                TomlFloat number => number.Value,
                TomlInt number => number.Value,
                _ => throw new MetaMorpheusException($"{key} must be a number, not {value.ReadableTypeName} '{value}'.")
            };
            if (key == LegacyQValueThresholdKey)
            {
                _legacyQValueThreshold = threshold;
            }
            else
            {
                _legacyPepQValueThreshold = threshold;
            }
            _oldTaskFile = true;
            return true;
        }

        /// <summary>The QValueThreshold the task TOML this was read from held; null when it held none.</summary>
        [TomlIgnore] public double? LegacyQValueThreshold => _legacyQValueThreshold;

        /// <summary>The PepQValueThreshold the task TOML this was read from held; null when it held none.</summary>
        [TomlIgnore] public double? LegacyPepQValueThreshold => _legacyPepQValueThreshold;

        /// <summary>
        /// True when <see cref="IdentificationFilter"/> was read from a task TOML that names neither the filter nor
        /// <see cref="FilterThreshold"/> and holds the old PEP thresholds (PepQValueThreshold below QValueThreshold, as the
        /// old GUI's PEP box saved them): the filter is then the PEP q-value at that PEP threshold, and the search says so
        /// once.
        /// </summary>
        [TomlIgnore]
        public bool IdentificationFilterFromLegacyThresholds =>
            !_identificationFilterNamed && !_filterThresholdNamed && _oldTaskFile
            && (_legacyPepQValueThreshold ?? LegacyPepQValueThresholdDefault) < (_legacyQValueThreshold ?? LegacyQValueThresholdDefault);

        // A copy keeps whether the filter and its threshold were named, or are read from old thresholds, so it resolves
        // them as the original does and the search still says where they came from.
        private void CopyIdentificationFilterRecord(CommonParameters source)
        {
            _identificationFilter = source._identificationFilter;
            _identificationFilterNamed = source._identificationFilterNamed;
            _filterThreshold = source._filterThreshold;
            _filterThresholdNamed = source._filterThresholdNamed;
            _legacyQValueThreshold = source._legacyQValueThreshold;
            _legacyPepQValueThreshold = source._legacyPepQValueThreshold;
            _oldTaskFile = source._oldTaskFile;
        }
        public double ScoreCutoff { get; private set; }
        /// <summary>
        /// This parameter determines which PSMs/Peptides will be used as postive training examples
        /// when training the GBDT model for PEP. 
        /// </summary>
        public double QValueCutoffForPepCalculation { get; set; }
        public IDigestionParams DigestionParams { get; private set; }
        public bool ReportAllAmbiguity { get; private set; }
        public int? NumberOfPeaksToKeepPerWindow { get; private set; }
        public double? MinimumAllowedIntensityRatioToBasePeak { get; private set; }
        public double? WindowWidthThomsons { get; private set; }
        public int? NumberOfWindows { get; private set; }
        public bool NormalizePeaksAccrossAllWindows { get; private set; }
        public bool TrimMs1Peaks { get; private set; }
        public bool TrimMsMsPeaks { get; private set; }
        public List<ProductType> CustomIons { get; private set; }
        public bool AssumeOrphanPeaksAreZ1Fragments { get; private set; }
        public int MaxHeterozygousVariants { get; private set; }
        public int MinVariantDepth { get; private set; }
        public bool AddTruncations { get; private set; }
        public DissociationType DissociationType { get; private set; }
        public string SeparationType { get; private set; }
        public DoubleRange RetentionTimeRange { get; private set; }

        public DissociationType MS2ChildScanDissociationType { get; set; }
        public DissociationType MS3ChildScanDissociationType { get; set; }

        public bool UseMostAbundantPrecursorIntensity { get; set; }

        /// <summary>
        /// Which precursor mass is used to select theoretical proteoform candidates during search.
        /// Defaults to <see cref="EngineLayer.PrecursorMassMatchMode.Monoisotopic"/>.
        /// </summary>
        public PrecursorMassMatchMode PrecursorMassMatchMode { get; set; }
        public DIAparameters? DIAparameters { get; set; } //only for DIA analysis involving pseudo ms2 scan generation
        public IFragmentationParams FragmentationParameters { get; set; }
        public string RTPredictorName { get; private set; }

        /// <summary>
        /// The 5' cap applied to every RNA loaded from a database (e.g. O-3P-1).
        /// Null in protein mode, and in RNA mode when the mzLib default should be used.
        /// </summary>
        public ChemicalFormula FivePrimeTerminus { get; private set; }

        /// <summary>
        /// The 3' cap applied to every RNA loaded from a database (e.g. OH).
        /// Null in protein mode, and in RNA mode when the mzLib default should be used.
        /// </summary>
        public ChemicalFormula ThreePrimeTerminus { get; private set; }

        /// <summary>
        /// Parses a chemical formula, returning null (the mzLib default cap) if the text is blank or cannot be understood.
        /// </summary>
        public static ChemicalFormula TryParseFormula(string formula)
        {
            if (string.IsNullOrWhiteSpace(formula))
                return null;
            try
            {
                // An empty formula ("C0", "H-1H") is the same as no cap. Test the text, not AtomCount: it is a signed
                // sum, so valid caps that net to zero atoms ("H-1Na", "O-1S") must be kept.
                var parsed = ChemicalFormula.ParseFormula(formula.Trim());
                return parsed.Formula.Length == 0 ? null : parsed;
            }
            catch (Exception)
            {
                return null;
            }
        }

        public CommonParameters Clone()
        {
            CommonParameters c = new CommonParameters();
            foreach (PropertyInfo property in typeof(CommonParameters).GetProperties())
            {
                if (property.CanWrite) property.SetValue(c, property.GetValue(this));
            }
            c.CopyIdentificationFilterRecord(this);
            return c;
        }

        public CommonParameters CloneWithNewDissociationType(DissociationType dissociationType)
        {
            CommonParameters c = new CommonParameters();
            foreach (PropertyInfo property in typeof(CommonParameters).GetProperties())
            {
                if (property.CanWrite) property.SetValue(c, property.GetValue(this));
            }
            c.CopyIdentificationFilterRecord(this);
            c.DissociationType = dissociationType;
            return c;
        }

        public CommonParameters CloneWithNewTerminus(FragmentationTerminus? terminus = null, bool? addCompIons = null) //for use with speedy semi-specific searches to get both termini
        {
            if (terminus == null)
            {
                terminus = DigestionParams.FragmentationTerminus;
            }
            if (addCompIons == null)
            {
                addCompIons = AddCompIons;
            }
            return new CommonParameters(
                                TaskDescriptor,
                                DissociationType,                   
                                MS2ChildScanDissociationType,
                                MS3ChildScanDissociationType,
                                SeparationType,
                                DoPrecursorDeconvolution,
                                UseProvidedPrecursorInfo,
                                DeconvolutionIntensityRatio,
                                DeconvolutionMaxAssumedChargeState,
                                ReportAllAmbiguity,
                                addCompIons.Value,//possibly changed
                                TotalPartitions,
                                IdentificationFilter,
                                FilterThreshold,
                                QValueCutoffForPepCalculation,
                                ScoreCutoff,
                                NumberOfPeaksToKeepPerWindow,
                                MinimumAllowedIntensityRatioToBasePeak,
                                WindowWidthThomsons,
                                NumberOfWindows,
                                NormalizePeaksAccrossAllWindows,
                                TrimMs1Peaks,
                                TrimMsMsPeaks,
                                ProductMassTolerance,
                                PrecursorMassTolerance,
                                ProductMassTolerance_LowRes,
                                DeconvolutionMassTolerance,
                                MaxThreadsToUsePerFile,
                                DigestionParams.Clone(terminus),
                                ListOfModsVariable,
                                ListOfModsFixed,
                                AssumeOrphanPeaksAreZ1Fragments,
                                MaxHeterozygousVariants,
                                MinVariantDepth,
                                AddTruncations,
                                PrecursorDeconvolutionParameters, 
                                ProductDeconvolutionParameters,
                                UseMostAbundantPrecursorIntensity,
                                DIAparameters,
                                FragmentationParameters,
                                PrecursorMassMatchMode,
                                RTPredictorName,
                                 RetentionTimeRange,
                                 FivePrimeTerminus?.Formula,
                                 ThreePrimeTerminus?.Formula);
        }

        /// <summary>
        /// Copy with a different TotalPartitions. Returns a new instance rather than mutating, because
        /// SetAllFileSpecificCommonParams hands back the task's own CommonParameters when a file has no
        /// file-specific settings — mutating that would rewrite the settings the task reports to the user.
        /// DigestionParams is shared rather than cloned so the digestion identity peptides carry is unchanged.
        /// </summary>
        public CommonParameters CloneWithNewTotalPartitions(int totalPartitions)
        {
            CommonParameters clone = new CommonParameters(
                                TaskDescriptor,
                                DissociationType,
                                MS2ChildScanDissociationType,
                                MS3ChildScanDissociationType,
                                SeparationType,
                                DoPrecursorDeconvolution,
                                UseProvidedPrecursorInfo,
                                DeconvolutionIntensityRatio,
                                DeconvolutionMaxAssumedChargeState,
                                ReportAllAmbiguity,
                                AddCompIons,
                                totalPartitions, //changed
                                IdentificationFilter,
                                FilterThreshold,
                                QValueCutoffForPepCalculation,
                                ScoreCutoff,
                                NumberOfPeaksToKeepPerWindow,
                                MinimumAllowedIntensityRatioToBasePeak,
                                WindowWidthThomsons,
                                NumberOfWindows,
                                NormalizePeaksAccrossAllWindows,
                                TrimMs1Peaks,
                                TrimMsMsPeaks,
                                ProductMassTolerance,
                                PrecursorMassTolerance,
                                ProductMassTolerance_LowRes,
                                DeconvolutionMassTolerance,
                                MaxThreadsToUsePerFile,
                                DigestionParams,
                                ListOfModsVariable,
                                ListOfModsFixed,
                                AssumeOrphanPeaksAreZ1Fragments,
                                MaxHeterozygousVariants,
                                MinVariantDepth,
                                AddTruncations,
                                PrecursorDeconvolutionParameters,
                                ProductDeconvolutionParameters,
                                UseMostAbundantPrecursorIntensity,
                                DIAparameters,
                                FragmentationParameters,
                                PrecursorMassMatchMode,
                                RTPredictorName,
                                 RetentionTimeRange,
                                 FivePrimeTerminus?.Formula,
                                 ThreePrimeTerminus?.Formula);

            // CustomIons is not a constructor parameter — the constructor reads it from the global
            // dissociation-type dictionary — so copy it across explicitly. GlycoSearchEngine branches on
            // CommonParameters.CustomIons, and this clone is handed to it.
            clone.CustomIons = CustomIons;
            return clone;
        }

        /// <summary>
        /// Copy with a different MaxThreadsToUsePerFile, for dividing a task's thread budget across spectra files searched at
        /// once. Returns a new instance for the same reason as <see cref="CloneWithNewTotalPartitions"/>: the task's own
        /// CommonParameters may be what a file was handed. Construct copies on one thread: the constructor reads and clears the
        /// process-wide custom product-type list.
        /// </summary>
        public CommonParameters CloneWithNewMaxThreadsToUsePerFile(int maxThreadsToUsePerFile)
        {
            // Copied by the method that already copies everything, then the one setting that differs is set on the copy - a
            // second copy of the constructor's argument list would be one more place to forget when CommonParameters gains a
            // setting. The setter is private, which is no obstacle from inside the class.
            CommonParameters clone = CloneWithNewTotalPartitions(TotalPartitions);
            clone.MaxThreadsToUsePerFile = maxThreadsToUsePerFile;
            return clone;
        }

        public void SetCustomProductTypes()
        {
            DigestionParams.ProductsFromDissociationType()[MassSpectrometry.DissociationType.Custom] = CustomIons;
        }

        public AnalyteType DetermineAnalyteType()
        {
            return DigestionParams switch
            {
                RnaDigestionParams => AnalyteType.Oligo,
                DigestionParams { Protease: not null } when DigestionParams.DigestionAgent.Name == "top-down"
                    => AnalyteType.Proteoform,
                _ => AnalyteType.Peptide
            };
        }
    }
}
