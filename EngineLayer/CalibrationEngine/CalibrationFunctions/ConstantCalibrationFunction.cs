using MathNet.Numerics.Statistics;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using System.Collections.Generic;

namespace EngineLayer.Calibration
{
    public class ConstantCalibrationFunction : ILearner<double>
    {
        #region Public Methods

        public IPredictorModel<double> Learn(F64Matrix observations, double[] targets)
        {
            double a = targets.Median();
            return new ConstantCalibrationFunctionPredictorModel(a);
        }

        #endregion Public Methods
    }

    internal class ConstantCalibrationFunctionPredictorModel : IPredictorModel<double>
    {
        #region Private Fields

        private double a;

        #endregion Private Fields

        #region Public Constructors

        public ConstantCalibrationFunctionPredictorModel(double a)
        {
            this.a = a;
        }

        #endregion Public Constructors

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
            return a;
        }

        #endregion Public Methods
    }
}