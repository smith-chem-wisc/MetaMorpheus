using MathNet.Numerics;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using System;

namespace EngineLayer.Calibration
{
    public class MzMultiplier : ILearner<double>
    {
        #region Public Methods

        public IPredictorModel<double> Learn(F64Matrix observations, double[] targets)
        {
            double[][] ok = new double[targets.Length][];
            for (int i = 0; i < targets.Length; i++)
            {
                ok[i] = observations.Row(i);
            }

            var ye = new Func<double[], double>[1];

            ye[0] = a => a[0];

            var f = Fit.LinearMultiDimFunc(ok, targets, ye);

            return new LinearCalibrationFunctionPredictorModel(f);
        }

        #endregion Public Methods
    }
}