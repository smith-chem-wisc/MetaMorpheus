using SharpLearning.Common.Interfaces;

namespace EngineLayer.Calibration
{
    internal class MS2Predictor
    {
        #region Private Fields

        private readonly IPredictorModel<double> bestMS2predictor;

        #endregion Private Fields

        #region Public Constructors

        public MS2Predictor(IPredictorModel<double> bestMS2predictor)
        {
            this.bestMS2predictor = bestMS2predictor;
        }

        #endregion Public Constructors

        #region Internal Methods

        internal double Predict(double mz, double retentionTime, double LOGtotalIonCurrent, double LOGinjectionTime, double logIntensity)
        {
            return bestMS2predictor.Predict(new[] { mz, retentionTime, LOGtotalIonCurrent, LOGinjectionTime, logIntensity });
        }

        #endregion Internal Methods
    }
}