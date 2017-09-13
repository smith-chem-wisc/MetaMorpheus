using SharpLearning.Common.Interfaces;
using System.Collections.Generic;

namespace EngineLayer.Calibration
{
    public class IdentityCalibrationFunctionPredictorModel : IPredictorModel<double>
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

        #endregion Public Methods
    }
}