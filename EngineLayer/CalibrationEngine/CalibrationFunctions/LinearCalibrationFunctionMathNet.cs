using MathNet.Numerics;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using System;
using System.Collections.Generic;

namespace EngineLayer.Calibration
{
    public class LinearCalibrationFunctionMathNet : ILearner<double>
    {
        #region Private Fields

        private readonly TransformFunction transformFunction;

        #endregion Private Fields

        #region Public Constructors

        public LinearCalibrationFunctionMathNet(TransformFunction transformFunction)
        {
            this.transformFunction = transformFunction;
        }

        public LinearCalibrationFunctionMathNet()
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public IPredictorModel<double> Learn(F64Matrix observations, double[] targets)
        {
            var numFeatures = transformFunction.numOutputs;
            double[][] ok = new double[targets.Length][];
            int k = 0;
            for (int i = 0; i < targets.Length; i++)
            {
                ok[k] = transformFunction.Transform(observations.Row(i));
                k++;
            }

            var ye = new Func<double[], double>[numFeatures + 1];
            ye[0] = a => 1;
            for (int i = 0; i < numFeatures; i++)
            {
                int j = i;
                ye[j + 1] = a => a[j];
            }
            var f = Fit.LinearMultiDimFunc(ok, targets, ye);

            return new LinearCalibrationFunctionPredictorModel(transformFunction, f);
        }

        #endregion Public Methods
    }

    internal class LinearCalibrationFunctionPredictorModel : IPredictorModel<double>
    {
        #region Private Fields

        private readonly TransformFunction transformFunction;
        private readonly Func<double[], double> f;

        #endregion Private Fields

        #region Public Constructors

        public LinearCalibrationFunctionPredictorModel(TransformFunction transformFunction, Func<double[], double> f)
        {
            this.transformFunction = transformFunction;
            this.f = f;
        }

        #endregion Public Constructors

        #region Public Methods

        public double[] GetRawVariableImportance()
        {
            throw new NotImplementedException();
        }

        public Dictionary<string, double> GetVariableImportance(Dictionary<string, int> featureNameToIndex)
        {
            throw new NotImplementedException();
        }

        public double Predict(double[] observation)
        {
            return f(transformFunction.Transform(observation));
        }

        #endregion Public Methods
    }
}