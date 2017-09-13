using SharpLearning.Common.Interfaces;

namespace EngineLayer.Calibration
{
    internal class MS1predictor
    {
        #region Private Fields

        private IPredictorModel<double> bestMS1predictor;

        #endregion Private Fields

        #region Public Constructors

        public MS1predictor(IPredictorModel<double> bestMS1predictor)
        {
            this.bestMS1predictor = bestMS1predictor;
        }

        #endregion Public Constructors

        #region Internal Methods

        internal double Predict(double mz, double retentionTime, double totalIonCurrent, double injectionTime)
        {
            return bestMS1predictor.Predict(new[] { mz, retentionTime, totalIonCurrent, injectionTime });
        }

        #endregion Internal Methods
    }
}