using MzLibUtil;
using System.Collections.Generic;

namespace EngineLayer
{
    public interface ICommonParameters
    {
        #region Public Properties

        int MaxParallelFilesToAnalyze { get; }
        int MaxThreadsToUsePerFile { get; }
        bool LocalizeAll { get; }
        IEnumerable<(string, string)> ListOfModsFixed { get; }
        IEnumerable<(string, string)> ListOfModsVariable { get; }
        IEnumerable<string> ListOfModTypesLocalize { get; }
        bool DoPrecursorDeconvolution { get; }
        bool UseProvidedPrecursorInfo { get; }
        double DeconvolutionIntensityRatio { get; }
        int DeconvolutionMaxAssumedChargeState { get; }
        Tolerance DeconvolutionMassTolerance { get; }
        int TotalPartitions { get; }
        bool BIons { get; }
        bool YIons { get; }
        bool ZdotIons { get; }
        bool CIons { get; }
        Tolerance ProductMassTolerance { get; }
        Tolerance PrecursorMassTolerance { get; }
        bool ConserveMemory { get; }
        double ScoreCutoff { get; }
        IDigestionParams DigestionParams { get; }
        bool ReportAllAmbiguity { get; }
        int TopNpeaks { get; }
        double MinRatio { get; }
        bool TrimMs1Peaks { get; }
        bool TrimMsMsPeaks { get; }
        string TaskDescriptor { get; }
        bool UseDeltaScore { get; }
        bool CalculateEValue { get; }

        #endregion Public Properties
    }
}