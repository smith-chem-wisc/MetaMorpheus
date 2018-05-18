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


        public readonly int MaxThreadsToUsePerFile;
        public readonly IEnumerable<(string, string)> ListOfModsFixed;
        public readonly IEnumerable<(string, string)> ListOfModsVariable;
        public readonly bool DoPrecursorDeconvolution;
        public readonly bool UseProvidedPrecursorInfo;
        public readonly double DeconvolutionIntensityRatio;
        public readonly int DeconvolutionMaxAssumedChargeState;
        public readonly Tolerance DeconvolutionMassTolerance;

        public readonly int TotalPartitions;

        public readonly bool BIons;

        public readonly bool YIons;

        public readonly bool ZdotIons;

        public readonly bool CIons;

        public readonly Tolerance ProductMassTolerance;
        public readonly Tolerance PrecursorMassTolerance;

        public readonly bool CompIons;

        public readonly double ScoreCutoff;

        public readonly DigestionParams DigestionParams;

        public readonly bool ReportAllAmbiguity;

        public readonly int TopNpeaks;
        public readonly double MinRatio;
        public readonly bool TrimMs1Peaks;
        public readonly bool TrimMsMsPeaks;
        public readonly string TaskDescriptor;

        public readonly bool UseDeltaScore;
        public readonly bool CalculateEValue;


        #region Public Constructors

        public CommonParameters(bool BIons = true, bool YIons = true, bool ZdotIons = false, bool CIons = false, bool DoPrecursorDeconvolution = true,
            bool UseProvidedPrecursorInfo = true, double DeconvolutionIntensityRatio = 3, int DeconvolutionMaxAssumedChargeState = 12, bool ReportAllAmbiguity = true,
            bool CompIons = false, int TotalPartitions = 1, double ScoreCutoff = 5, int TopNpeaks = 200, double MinRatio = 0.01, bool TrimMs1Peaks = false,
            bool TrimMsMsPeaks = true, bool UseDeltaScore = false, bool CalculateEValue = false, int prodMassTol = 20, int preMassTol = 5, int deconMassTol = 4,
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
            ProductMassTolerance = new PpmTolerance(prodMassTol); 
            PrecursorMassTolerance = new PpmTolerance(preMassTol);
            DeconvolutionMassTolerance = new PpmTolerance(deconMassTol);

            if (DigestionParams == null)
            {
                this.DigestionParams = new DigestionParams(); // not compile time constant
            }
            if (ListOfModsVariable ==null)
            {
                this.ListOfModsVariable = new List<(string, string)> { ("Common Variable", "Oxidation of M") };
            }
            if (ListOfModsFixed == null)
            {
                this.ListOfModsFixed = new List<(string, string)> { ("Common Fixed", "Carbamidomethyl of C"), ("Common Fixed", "Carbamidomethyl of U") };
            }
            

            
            
            
           
        }

        #endregion Public Constructors


        #region Public Methods

        public CommonParameters Clone()
        {
            return (CommonParameters)this.MemberwiseClone();
        }

        public ParallelOptions ParallelOptions()
        {
            return new ParallelOptions { MaxDegreeOfParallelism = MaxParallelFilesToAnalyze };
        }

        #endregion Public Methods
    }
}