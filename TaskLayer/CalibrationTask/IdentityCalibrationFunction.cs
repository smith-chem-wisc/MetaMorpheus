using EngineLayer.Calibration;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;

namespace TaskLayer
{
    internal class IdentityCalibrationFunction : ILearner<double>
    {
        #region Public Methods

        public IPredictorModel<double> Learn(F64Matrix observations, double[] targets)
        {
            return new IdentityCalibrationFunctionPredictorModel();
        }

        #endregion Public Methods
    }
}