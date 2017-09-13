using SharpLearning.Common.Interfaces;

namespace EngineLayer.Calibration
{
    internal class MS2predictor
    {
        #region Private Fields

        private IPredictorModel<double> bestMS2predictor;

        #endregion Private Fields

        #region Public Constructors

        public MS2predictor(IPredictorModel<double> bestMS2predictor)
        {
            this.bestMS2predictor = bestMS2predictor;
        }

        #endregion Public Constructors

        #region Internal Methods

        internal double Predict(double mz, double retentionTime, double totalIonCurrent, double injectionTime)
        {
            return bestMS2predictor.Predict(new[] { mz, retentionTime, totalIonCurrent, injectionTime });
        }

        #endregion Internal Methods
    }
}