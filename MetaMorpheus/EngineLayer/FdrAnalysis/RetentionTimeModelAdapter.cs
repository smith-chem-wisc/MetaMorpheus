using PredictionClients.Koina.AbstractClasses;
using System.Collections.Generic;

namespace EngineLayer.FdrAnalysis
{
    /// <summary>
    /// Wraps a mzLib RetentionTimeModel (from the NuGet package) to implement the
    /// IRetentionTimePredictor interface used by PepAnalysisEngine.
    /// This adapter exists because we cannot modify the mzLib NuGet source to add the interface directly.
    /// </summary>
    public class RetentionTimeModelAdapter : IRetentionTimePredictor
    {
        private readonly RetentionTimeModel _model;

        public RetentionTimeModelAdapter(RetentionTimeModel model)
        {
            _model = model;
        }

        public List<PeptideRTPrediction> Predict(List<RetentionTimePredictionInput> inputs)
            => _model.Predict(inputs);

        public bool IsIndexedRetentionTimeModel => _model.IsIndexedRetentionTimeModel;

        public string ModelName => _model.ModelName;
    }
}
