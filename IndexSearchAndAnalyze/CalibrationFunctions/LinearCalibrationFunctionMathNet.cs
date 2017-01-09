using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace IndexSearchAndAnalyze
{
    public class LinearCalibrationFunctionMathNet : CalibrationFunction
    {
        private Func<double[], double> f;
        private int numFeatures;
        private TransformFunction transformFunction;

        public LinearCalibrationFunctionMathNet(TransformFunction transformFunction)
        {
            this.transformFunction = transformFunction;
            numFeatures = transformFunction.numOutputs;
        }

        public override double Predict(double[] t)
        {
            return f(transformFunction.Transform(t));
        }

        public override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            double[][] ok = new double[trainingList.Count()][];
            int k = 0;
            foreach (LabeledDataPoint p in trainingList)
            {
                ok[k] = transformFunction.Transform(p.inputs);
                k++;
            }
            var ok2 = trainingList.Select(b => b.output).ToArray();

            var ye = new Func<double[], double>[numFeatures + 1];
            ye[0] = a => 1;
            for (int i = 0; i < numFeatures; i++)
            {
                int j = i;
                ye[j + 1] = a => a[j];
            }
            f = Fit.LinearMultiDimFunc(ok, ok2, ye);
            //onOutput("Finished fitting a linear"));
        }
    }
}