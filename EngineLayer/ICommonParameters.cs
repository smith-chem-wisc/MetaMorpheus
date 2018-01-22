using MzLibUtil;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public interface ICommonParameters
    {
        #region Public Properties

        int? MaxParallelFilesToAnalyze { get; }
        int MaxThreadsToUsePerFile { get; }
        bool LocalizeAll { get; }
        List<Tuple<string, string>> ListOfModsFixed { get; }
        List<Tuple<string, string>> ListOfModsVariable { get; }
        List<Tuple<string, string>> ListOfModsLocalize { get; }

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

        DigestionParams DigestionParams { get; }

        bool ReportAllAmbiguity { get; }

        bool ExcelCompatible { get; }
        int? TopNpeaks { get; }
        double? MinRatio { get; }
        bool TrimMs1Peaks { get; }
        bool TrimMsMsPeaks { get; }
        string TaskDescriptor { get; }

        bool CalculateEValue { get; }

        #endregion Public Properties
    }
}