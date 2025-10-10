using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Reflection;
using Nett;
using Omics.Digestion;
using Omics.Fragmentation.Peptide;
using Transcriptomics.Digestion;
using EngineLayer.DIA;

namespace EngineLayer
{
    public class CommonParameters
    {
        // this parameterless constructor needs to exist to read the toml.
        // if you can figure out a way to get rid of it, feel free...
        public CommonParameters()
            : this(digestionParams: null)
        {
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
            double qValueThreshold = 0.01,
            double pepQValueThreshold = 1.0,
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
            Tolerance deconvolutionMassTolerance = null,
            int maxThreadsToUsePerFile = -1, 
            IDigestionParams digestionParams = null, 
            IEnumerable<(string, string)> listOfModsVariable = null, 
            IEnumerable<(string, string)> listOfModsFixed = null,
            bool assumeOrphanPeaksAreZ1Fragments = true,
            int maxSequenceVariantsPerIsoform = 0,
            int minAlleleDepth = 0,
            int maxSequenceVariantIsoforms = 1,
            bool addTruncations = false,
            DeconvolutionParameters precursorDeconParams = null,
            DeconvolutionParameters productDeconParams = null,
            bool useMostAbundantPrecursorIntensity = true,
            DIAparameters diaParameters = null)

        {
            TaskDescriptor = taskDescriptor;
            DoPrecursorDeconvolution = doPrecursorDeconvolution;
            UseProvidedPrecursorInfo = useProvidedPrecursorInfo;
            DeconvolutionIntensityRatio = deconvolutionIntensityRatio;
            ReportAllAmbiguity = reportAllAmbiguity;
            AddCompIons = addCompIons;
            TotalPartitions = totalPartitions;
            QValueThreshold = qValueThreshold;
            PepQValueThreshold = pepQValueThreshold;
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
            DeconvolutionMassTolerance = deconvolutionMassTolerance ?? new PpmTolerance(4);
            DigestionParams = digestionParams ?? new DigestionParams();
            DissociationType = dissociationType;
            SeparationType = separationType;
            MS2ChildScanDissociationType = ms2childScanDissociationType;
            MS3ChildScanDissociationType = ms3childScanDissociationType;
            UseMostAbundantPrecursorIntensity = useMostAbundantPrecursorIntensity;
            AssumeOrphanPeaksAreZ1Fragments = assumeOrphanPeaksAreZ1Fragments;
            MaxSequenceVariantsPerIsoform = maxSequenceVariantsPerIsoform;
            MinAlleleDepth = minAlleleDepth;
            MaxSequenceVariantIsoforms = maxSequenceVariantIsoforms;
            AddTruncations = addTruncations;
            DIAparameters = diaParameters;

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
            }
            else
            {
                ListOfModsVariable = listOfModsVariable ?? new List<(string, string)> { ("Common Variable", "Oxidation on M") };
                ListOfModsFixed = listOfModsFixed ?? new List<(string, string)> { ("Common Fixed", "Carbamidomethyl on C"), ("Common Fixed", "Carbamidomethyl on U") };
            }

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
        public int MaxThreadsToUsePerFile { get; private set; }
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
        public Tolerance PrecursorMassTolerance { get; set; } // public setter required for calibration task
        public bool AddCompIons { get; private set; }
        /// <summary>
        /// Only peptides/PSMs with Q-Value and Q-Value Notch below this threshold are used for quantification and
        /// spectral library generation. If SearchParameters.WriteHighQValuePsms is set to false, only 
        /// peptides/PSMs with Q-Value/Notch below this threshold are included in .psmtsv results files.
        /// </summary>
        public double QValueThreshold { get; private set; }
        /// <summary>
        /// Only peptides/PSMs with PEP Q-Value below this threshold are used for quantification and
        /// spectral library generation. If SearchParameters.WriteHighQValuePsms is set to false, only 
        /// peptides/PSMs with PEP Q-Value below this threshold are included in .psmtsv results files.
        /// </summary>
        public double PepQValueThreshold { get; private set; }
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
        public int MaxSequenceVariantsPerIsoform { get; private set; }
        public int MinAlleleDepth { get; private set; }
        public int MaxSequenceVariantIsoforms { get; private set; }
        public bool AddTruncations { get; private set; }
        public DissociationType DissociationType { get; private set; }
        public string SeparationType { get; private set; }

        public DissociationType MS2ChildScanDissociationType { get; private set; }
        public DissociationType MS3ChildScanDissociationType { get; private set; }

        public bool UseMostAbundantPrecursorIntensity { get; set; }
        public DIAparameters? DIAparameters { get; set; } //only for DIA analysis involving pseudo ms2 scan generation

        public CommonParameters Clone()
        {
            CommonParameters c = new CommonParameters();
            foreach (PropertyInfo property in typeof(CommonParameters).GetProperties())
            {
                property.SetValue(c, property.GetValue(this));
            }
            return c;
        }

        public CommonParameters CloneWithNewDissociationType(DissociationType dissociationType)
        {
            CommonParameters c = new CommonParameters();
            foreach (PropertyInfo property in typeof(CommonParameters).GetProperties())
            {
                property.SetValue(c, property.GetValue(this));
            }
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
                                QValueThreshold,
                                PepQValueThreshold,
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
                                DeconvolutionMassTolerance,
                                MaxThreadsToUsePerFile,
                                DigestionParams.Clone(terminus),
                                ListOfModsVariable,
                                ListOfModsFixed,
                                AssumeOrphanPeaksAreZ1Fragments,
                                MaxSequenceVariantsPerIsoform,
                                MinAlleleDepth,
                                MaxSequenceVariantIsoforms,
                                AddTruncations,
                                PrecursorDeconvolutionParameters, 
                                ProductDeconvolutionParameters);
        }

        public void SetCustomProductTypes()
        {
            DigestionParams.ProductsFromDissociationType()[MassSpectrometry.DissociationType.Custom] = CustomIons;
        }
    }
}
