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
        public readonly bool UseDeltaScore;
        public readonly bool CalculateEValue;

        #region Public Constructors

        public CommonParameters(bool BIons = true, bool YIons = true, bool ZdotIons = false, bool CIons = false, bool DoPrecursorDeconvolution = true,
            bool UseProvidedPrecursorInfo = true, double DeconvolutionIntensityRatio = 3, int DeconvolutionMaxAssumedChargeState = 12, bool ReportAllAmbiguity = true,
            bool CompIons = false, int TotalPartitions = 1, double ScoreCutoff = 5, int TopNpeaks = 200, double MinRatio = 0.01, bool TrimMs1Peaks = false,
            bool TrimMsMsPeaks = true, bool UseDeltaScore = false, bool CalculateEValue = false, double prodMassTol = 20, double preMassTol = 5, double deconMassTol = 4,
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

        public CommonParameters Clone()
        {
            // same settings as this.[myparameters]
            return new CommonParameters();
        }

        public override string ToString()
        {
            string ret = "";

            ret += "BIons =" + BIons + ";";
            ret += "YIons =" + YIons + ";";
            ret += "ZdotIons =" + ZdotIons + ";";
            ret += "CIons =" + CIons + ";";
            ret += "DoPreCursorDeconvolution =" + DoPrecursorDeconvolution + ";";
            ret += "UseProvidedPrecursorInfo =" + UseProvidedPrecursorInfo + ";";
            ret += "DeconvolutionIntensityRatio =" + DeconvolutionIntensityRatio + ";";
            ret += "DeconvolutionMaxAssumedChargeState =" + DeconvolutionMaxAssumedChargeState + ";";
            ret += "ReportAllAmbiguity =" + ReportAllAmbiguity + ";";
            ret += "CompIons =" + CompIons + ";";
            ret += "TotalPartitions =" + TotalPartitions + ";";
            ret += "ScoreCutoff =" + ScoreCutoff + ";";
            ret += "TopNpeaks =" + TopNpeaks + ";";
            ret += "MinRatio =" + MinRatio + ";";
            ret += "TrimMs1Peaks =" + TrimMs1Peaks + ";";
            ret += "TrimMsMsPeaks =" + TrimMsMsPeaks + ";";
            ret += "UseDeltaScore =" + UseDeltaScore + ";";
            ret += "CalculateEValue =" + CalculateEValue + ";";
            ret += "prodMassTol =" + ProductMassTolerance.Value + ";";
            ret += "preMassTol =" + PrecursorMassTolerance.Value + ";";
            ret += "deconMassTol =" + DeconvolutionMassTolerance.Value + ";";
            ret += "MaxThreadsToUsePerFile =" + MaxThreadsToUsePerFile + ";";
            ret += "DigestionParams: Protease =" + DigestionParams.Protease.Name + ";";
            ret += "DigestionParams: MaxMissedCleavages =" + DigestionParams.MaxMissedCleavages + ";";
            ret += "DigestionParams: MinPeptideLength =" + DigestionParams.MinPeptideLength + ";";
            ret += "DigestionParams: MaxPeptideLength =" + DigestionParams.MaxPeptideLength + ";";
            ret += "DigestionParams: MaxModificationIsoforms =" + DigestionParams.MaxModificationIsoforms + ";";
            ret += "DigestionParams: InitiatorMethionineBehavior =" + DigestionParams.InitiatorMethionineBehavior + ";";
            ret += "DigestionParams: MaxModsForPeptides =" + DigestionParams.MaxModsForPeptide + ";";
            ret += "DigestionParams: SemiProteaseDigestion =" + DigestionParams.SemiProteaseDigestion + ";";
            ret += "DigestionParams: TerminusTypeSemiProtease =" + DigestionParams.TerminusTypeSemiProtease + ";";

            string contentVariable = "";
            foreach (var thing in ListOfModsVariable)
            {
               
                contentVariable+=  thing.Item1 + "/";
                contentVariable += thing.Item2 + ",";
            }
            ret += "ListOfModsVariable =" + contentVariable + ";";
            string contentFixed = "";
            foreach (var thing in ListOfModsFixed)
            {

                contentFixed += thing.Item1 + "/";
                contentFixed += thing.Item2 + ",";
            }
            ret += "ListOfModsFixed =" + contentFixed + ";";
     
            return ret;
        }

        #endregion Public Methods
    }
}
