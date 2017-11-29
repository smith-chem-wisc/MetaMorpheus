using SharpLearning.Common.Interfaces;

namespace EngineLayer.Calibration
{
    internal class MS1Predictor
    {
        #region Private Fields

        private readonly IPredictorModel<double> bestMS1predictor;

        #endregion Private Fields

        #region Public Constructors

        public MS1Predictor(IPredictorModel<double> bestMS1predictor)
        {
            this.bestMS1predictor = bestMS1predictor;
        }

        #endregion Public Constructors

        #region Internal Methods

        internal double Predict(double mz, double retentionTime, double logTotalIonCurrent, double logInjectionTime, double logIntensity)
        {
            return bestMS1predictor.Predict(new[] { mz, retentionTime, logTotalIonCurrent, logInjectionTime, logIntensity });
        }

        #endregion Internal Methods
    }
}