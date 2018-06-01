using EngineLayer;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class CommonParameters
    {
        //Any new property must also be added in MetaMorpheusTask.SetAllFileSpecificCommonParams, else it be overwritten by file specific params
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
            int MaxThreadsToUsePerFile = -1, DigestionParams DigestionParams = null, List<(string, string)> ListOfModsVariable = null, List<(string, string)> ListOfModsFixed = null)
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

            if (ProductMassTolerance == null)
            {
                this.ProductMassTolerance = new PpmTolerance(20);
            }
            else
            {
                this.ProductMassTolerance = ProductMassTolerance;
            }

            if (PrecursorMassTolerance == null)
            {
                this.PrecursorMassTolerance = new PpmTolerance(5);
            }
            else
            {
                this.PrecursorMassTolerance = PrecursorMassTolerance;
            }

            if (DeconvolutionMassTolerance == null)
            {
                this.DeconvolutionMassTolerance = new PpmTolerance(4);
            }
            else
            {
                this.DeconvolutionMassTolerance = DeconvolutionMassTolerance;
            }
            
            if (DigestionParams == null)
            {
                this.DigestionParams = new DigestionParams(); // not compile time constant
            }
            else
            {
                this.DigestionParams = DigestionParams;
            }

            if (ListOfModsVariable == null)
            {
                this.ListOfModsVariable = new List<(string, string)> { ("Common Variable", "Oxidation of M") };
            }
            else
            {
                this.ListOfModsVariable = ListOfModsVariable;
            }

            if (ListOfModsFixed == null)
            {
                this.ListOfModsFixed = new List<(string, string)> { ("Common Fixed", "Carbamidomethyl of C"), ("Common Fixed", "Carbamidomethyl of U") };
            }
            else
            {
                this.ListOfModsFixed = ListOfModsFixed;
            }
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

        #region Public Methods

        //public CommonParameters Clone()
        //{
        //    // same settings as this.[myparameters]
        //    return new CommonParameters();
        //}

        public void SetTaskDescriptor(string TaskDescriptor)
        {
            this.TaskDescriptor = TaskDescriptor;
        }

        private void SetMaxThreadsToUsePerFile(int MaxThreadsToUsePerFile)
        {
            this.MaxThreadsToUsePerFile = MaxThreadsToUsePerFile;
        }

        private void SetListOfModsFixed(IEnumerable<(string, string)> ListOfModsFixed)
        {
            this.ListOfModsFixed = ListOfModsFixed;
        }

        private void  SetListOfModsVariable(IEnumerable<(string, string)> ListOfModsVariable)
        {
            this.ListOfModsVariable = ListOfModsVariable;
        }

        private void SetDoPrecursorDeconvolution(bool DoPrecursorDeconvolution)
        {
            this.DoPrecursorDeconvolution = DoPrecursorDeconvolution;
        }

        private void SetUseProvidedPrecursorInfo(bool UseProvidedPrecursorInfo)
        {
            this.UseProvidedPrecursorInfo = UseProvidedPrecursorInfo;
        }

        private void SetDeconvolutionIntensityRatio(int DeconvolutionIntesityRatio)
        {
            this.DeconvolutionIntensityRatio = DeconvolutionIntensityRatio;
        }

        private void SetDeconvolutionMaxAssumedChargeState(int DeconvolutionMaxAssumedChargeState)
        {
            this.DeconvolutionMaxAssumedChargeState = DeconvolutionMaxAssumedChargeState;
        }

        public void SetDeconvolutionMassTolerance(Tolerance DeconvolutionMassTolerance)
        {
            this.DeconvolutionMassTolerance = DeconvolutionMassTolerance;
        }

        public void SetTotalPartitions(int TotalPartitions)
        {
            this.TotalPartitions = TotalPartitions;
        }

        public void SetBIons(bool BIons)
        {
            this.BIons = BIons;
        }

        public void SetYIons(bool YIons)
        {
            this.YIons = YIons;
        }

        public void SetZdotIons(bool ZdotIons)
        {
            this.ZdotIons = ZdotIons;
        }

        public void SetCIons(bool CIons)
        {
            this.CIons = CIons; 
        }

        public void SetProductMassTolerance(Tolerance ProductMassTolerance)
        {
            this.ProductMassTolerance = ProductMassTolerance;
        }

        public void SetPrecursorMassTolerance(Tolerance PrecursorMassTolerance)
        {
            this.PrecursorMassTolerance = PrecursorMassTolerance;
        }

        public void SetCompIons(bool CompIons)
        {
            this.CompIons = CompIons;
        }

        public void SetScoreCutoff(double ScoreCutoff)
        {
            this.ScoreCutoff = ScoreCutoff;
        }

        public void SetDigestionParams(DigestionParams DigetionsParams)
        {
            this.DigestionParams = DigestionParams;
        }

        public void SetReportAllAmbiguity(bool ReportAllAmbiguity)
        {
            this.ReportAllAmbiguity = ReportAllAmbiguity;
        }

        public void SetTopNpeaks(int TopNpeaks)
        {
            this.TopNpeaks = TopNpeaks;
        }

        public void SetMinRatio(double MinRatio)
        {
            this.MinRatio = MinRatio;
        }

        public void SetTrimMs1Peaks(bool TrimMS1Peaks)
        {
            this.TrimMs1Peaks = TrimMS1Peaks;
        }

        public void SetTrimMsMsPeaks(bool TrimMsMsPeaks)
        {
            this.TrimMsMsPeaks = TrimMsMsPeaks;
        }
        public void SetUseDeltaScore(bool UseDeltaScore)
        {
            this.UseDeltaScore = UseDeltaScore;
        }

        public void SetCalculateEValue(bool CalculateEValue)
        {
            this.CalculateEValue = CalculateEValue;
        }

        #endregion Public Methods
    }
}
