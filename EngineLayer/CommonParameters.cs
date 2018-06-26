using EngineLayer;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class CommonParameters
    {
        #region Public Constructors

        // this parameterless constructor needs to exist to read the toml. 
        // if you can figure out a way to get rid of it, feel free...
        public CommonParameters() : this(DigestionParams: null)
        {

        }

        public CommonParameters(bool BIons = true, bool YIons = true, bool ZdotIons = false, bool CIons = false, bool DoPrecursorDeconvolution = true,
            bool UseProvidedPrecursorInfo = true, double DeconvolutionIntensityRatio = 3, int DeconvolutionMaxAssumedChargeState = 12, bool ReportAllAmbiguity = true,
            bool CompIons = false, int TotalPartitions = 1, double ScoreCutoff = 5, int TopNpeaks = 200, double MinRatio = 0.01, bool TrimMs1Peaks = false,
            bool TrimMsMsPeaks = true, bool UseDeltaScore = false, bool CalculateEValue = false, Tolerance ProductMassTolerance = null, Tolerance PrecursorMassTolerance = null, Tolerance DeconvolutionMassTolerance = null,
            int MaxThreadsToUsePerFile = -1, DigestionParams DigestionParams = null, IEnumerable<(string, string)> ListOfModsVariable = null, IEnumerable<(string, string)> ListOfModsFixed = null)
        {
            this.BIons = BIons;
            this.YIons = YIons;
            this.ZdotIons = ZdotIons;
            this.CIons = CIons;
            this.DoPrecursorDeconvolution = DoPrecursorDeconvolution;
            this.UseProvidedPrecursorInfo = UseProvidedPrecursorInfo;
            this.DeconvolutionIntensityRatio = DeconvolutionIntensityRatio;
            this.DeconvolutionMaxAssumedChargeState = DeconvolutionMaxAssumedChargeState;
            this.ReportAllAmbiguity = ReportAllAmbiguity;
            this.CompIons = CompIons;
            this.TotalPartitions = TotalPartitions;
            this.ScoreCutoff = ScoreCutoff;
            this.TopNpeaks = TopNpeaks;
            this.MinRatio = MinRatio;
            this.TrimMs1Peaks = TrimMs1Peaks;
            this.TrimMsMsPeaks = TrimMsMsPeaks;
            this.UseDeltaScore = UseDeltaScore;
            this.CalculateEValue = CalculateEValue;
            this.MaxThreadsToUsePerFile = MaxThreadsToUsePerFile;

            this.ProductMassTolerance = ProductMassTolerance ?? new PpmTolerance(20);
            this.PrecursorMassTolerance = PrecursorMassTolerance ?? new PpmTolerance(5);
            this.DeconvolutionMassTolerance = DeconvolutionMassTolerance ?? new PpmTolerance(4);
            this.DigestionParams = DigestionParams ?? new DigestionParams();
            this.ListOfModsVariable = ListOfModsVariable ?? new List<(string, string)> { ("Common Variable", "Oxidation of M") };
            this.ListOfModsFixed = ListOfModsFixed ?? new List<(string, string)> { ("Common Fixed", "Carbamidomethyl of C"), ("Common Fixed", "Carbamidomethyl of U") };

            if (MaxThreadsToUsePerFile == -1)
            {
                this.MaxThreadsToUsePerFile = Environment.ProcessorCount > 1 ? Environment.ProcessorCount - 1 : 1;
            }
            else
            {
                this.MaxThreadsToUsePerFile = MaxThreadsToUsePerFile;
            }
        }

        #endregion Public Constructors

        #region Public Properties

        //Any new property must not be nullable (int?) or else if it is null, the null setting will not be written to a toml and the default will override (so it's okay if the default is null)
        public string TaskDescriptor { get; set; }
        public int MaxThreadsToUsePerFile { get; private set; }
        public IEnumerable<(string, string)> ListOfModsFixed { get; private set; }
        public IEnumerable<(string, string)> ListOfModsVariable { get; private set; }
        public bool DoPrecursorDeconvolution { get; private set; }
        public bool UseProvidedPrecursorInfo { get; private set; }
        public double DeconvolutionIntensityRatio { get; private set; }
        public int DeconvolutionMaxAssumedChargeState { get; private set; }
        public Tolerance DeconvolutionMassTolerance { get; private set; }
        public int TotalPartitions { get; private set; }
        public bool BIons { get; private set; }
        public bool YIons { get; private set; }
        public bool ZdotIons { get; private set; }
        public bool CIons { get; private set; }
        public Tolerance ProductMassTolerance { get; private set; }
        public Tolerance PrecursorMassTolerance { get; private set; }
        public bool CompIons { get; private set; }
        public double ScoreCutoff { get; private set; }
        public DigestionParams DigestionParams { get; private set; }
        public bool ReportAllAmbiguity { get; private set; }
        public int TopNpeaks { get; private set; }
        public double MinRatio { get; private set; }
        public bool TrimMs1Peaks { get; private set; }
        public bool TrimMsMsPeaks { get; private set; }
        public bool UseDeltaScore { get; private set; }
        public bool CalculateEValue { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public CommonParameters Clone()
        {
            return new CommonParameters(
                BIons: this.BIons,
                YIons: this.YIons,
                ZdotIons: this.ZdotIons,
                CIons: this.CIons,
                DoPrecursorDeconvolution: this.DoPrecursorDeconvolution,
                UseProvidedPrecursorInfo: this.UseProvidedPrecursorInfo,
                DeconvolutionIntensityRatio: this.DeconvolutionIntensityRatio,
                DeconvolutionMaxAssumedChargeState: this.DeconvolutionMaxAssumedChargeState,
                ReportAllAmbiguity: this.ReportAllAmbiguity,
                CompIons: this.CompIons,
                TotalPartitions: this.TotalPartitions,
                ScoreCutoff: this.ScoreCutoff,
                TopNpeaks: this.TopNpeaks,
                MinRatio: this.MinRatio,
                TrimMs1Peaks: this.TrimMs1Peaks,
                TrimMsMsPeaks: this.TrimMsMsPeaks,
                UseDeltaScore: this.UseDeltaScore,
                CalculateEValue: this.CalculateEValue,
                ProductMassTolerance: this.ProductMassTolerance,
                PrecursorMassTolerance: this.PrecursorMassTolerance,
                DeconvolutionMassTolerance: this.DeconvolutionMassTolerance,
                MaxThreadsToUsePerFile: this.MaxThreadsToUsePerFile,
                DigestionParams: this.DigestionParams,
                ListOfModsVariable: this.ListOfModsVariable,
                ListOfModsFixed: this.ListOfModsFixed
            );
        }
        
        public void SetProductMassTolerance(Tolerance ProductMassTolerance)
        {
            this.ProductMassTolerance = ProductMassTolerance;
        }

        public void SetPrecursorMassTolerance(Tolerance PrecursorMassTolerance)
        {
            this.PrecursorMassTolerance = PrecursorMassTolerance;
        }

        public void SetDigestionParams(DigestionParams DigestionParams)
        {
            this.DigestionParams = DigestionParams;
        }

        #endregion Public Methods
    }
}
