using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using System.Collections.Generic;

namespace EngineLayer.Calibration
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

    internal class IdentityCalibrationFunctionPredictorModel : IPredictorModel<double>
    {
        #region Public Methods

        public double[] GetRawVariableImportance()
        {
            throw new System.NotImplementedException();
        }

        public Dictionary<string, double> GetVariableImportance(Dictionary<string, int> featureNameToIndex)
        {
            throw new System.NotImplementedException();
        }

        public double Predict(double[] observation)
        {
            return 0;
        }

        public double[] Predict(F64Matrix observations)
        {
            throw new System.NotImplementedException();
        }

        #endregion Public Methods
    }
}