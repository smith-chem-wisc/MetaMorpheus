using MathNet.Numerics;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Calibration
{
    public class LinearCalibrationFunctionMathNet : ILearner<double>
    {
        #region Private Fields

        private readonly int[] v;

        #endregion Private Fields

        #region Public Constructors

        public LinearCalibrationFunctionMathNet(int[] v)
        {
            this.v = v;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return "LinearCalibrationFunctionMathNet " + string.Join(",", v);
        }

        public IPredictorModel<double> Learn(F64Matrix observations, double[] targets)
        {
            double[][] ok = new double[targets.Length][];
            for (int i = 0; i < targets.Length; i++)
            {
                ok[i] = observations.Row(i);
            }

            var ye = new Func<double[], double>[v.Length + 1];
            int k = 0;
            ye[0] = a => 1;
            for (int i = 0; i < observations.Row(0).Length; i++)
            {
                if (v.Contains(i))
                {
                    int ii = i;
                    int kk = k;
                    ye[kk + 1] = a => a[ii];
                    k++;
                }
            }
            var f = Fit.LinearMultiDimFunc(ok, targets, ye);

            return new LinearCalibrationFunctionPredictorModel(f);
        }

        #endregion Public Methods
    }

    internal class LinearCalibrationFunctionPredictorModel : IPredictorModel<double>
    {
        #region Private Fields

        private readonly Func<double[], double> f;

        #endregion Private Fields

        #region Public Constructors

        public LinearCalibrationFunctionPredictorModel(Func<double[], double> f)
        {
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
            return new Dictionary<string, double>();
        }

        public double Predict(double[] observation)
        {
            return f(observation);
        }

        public double[] Predict(F64Matrix observations)
        {
            throw new NotImplementedException();
        }

        #endregion Public Methods
    }
}