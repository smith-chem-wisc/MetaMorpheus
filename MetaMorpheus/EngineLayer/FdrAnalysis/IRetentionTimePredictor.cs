using PredictionClients.Koina.AbstractClasses;
using System.Collections.Generic;

namespace EngineLayer.FdrAnalysis
{
    /// <summary>
    /// Abstraction over any retention time prediction model.
    /// Allows PepAnalysisEngine to accept any RT predictor and enables
    /// mock injection in unit tests without requiring HTTP calls.
    /// </summary>
    public interface IRetentionTimePredictor
    {
        /// <summary>
        /// Predict retention time values for a list of peptide sequences.
        /// Input sequences should be in mzLib full-sequence format.
        /// Implementations must handle batching internally.
        /// </summary>
        List<PeptideRTPrediction> Predict(List<RetentionTimePredictionInput> inputs);

        /// <summary>
        /// True if this model predicts indexed retention time (iRT, dimensionless scale).
        /// False if the model predicts absolute retention time in minutes.
        /// When true, per-file linear calibration (iRT → RT) is applied.
        /// </summary>
        bool IsIndexedRetentionTimeModel { get; }

        /// <summary>
        /// Human-readable model name used in log messages and warnings.
        /// </summary>
        string ModelName { get; }
    }
}
