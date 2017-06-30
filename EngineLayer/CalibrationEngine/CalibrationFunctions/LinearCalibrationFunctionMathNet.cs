using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Calibration
{
    internal class LinearCalibrationFunctionMathNet : CalibrationFunction
    {
        #region Private Fields

        private readonly int numFeatures;
        private readonly TransformFunction transformFunction;
        private Func<double[], double> f;

        #endregion Private Fields

        #region Public Constructors

        public LinearCalibrationFunctionMathNet(TransformFunction transformFunction)
        {
            this.transformFunction = transformFunction;
            numFeatures = transformFunction.numOutputs;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("Linear");
            sb.Append(" numFeatures: " + numFeatures);
            sb.Append(" transform num outputs: " + transformFunction.numOutputs);
            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal override double Predict(double[] t)
        {
            return f(transformFunction.Transform(t));
        }

        internal override void Train<LabeledDataPoint>(IEnumerable<LabeledDataPoint> trainingList)
        {
            double[][] ok = new double[trainingList.Count()][];
            int k = 0;
            foreach (LabeledDataPoint p in trainingList)
            {
                ok[k] = transformFunction.Transform(p.Inputs);
                k++;
            }
            var ok2 = trainingList.Select(b => b.Label).ToArray();

            var ye = new Func<double[], double>[numFeatures + 1];
            ye[0] = a => 1;
            for (int i = 0; i < numFeatures; i++)
            {
                int j = i;
                ye[j + 1] = a => a[j];
            }
            f = Fit.LinearMultiDimFunc(ok, ok2, ye);
        }

        #endregion Internal Methods
    }
}