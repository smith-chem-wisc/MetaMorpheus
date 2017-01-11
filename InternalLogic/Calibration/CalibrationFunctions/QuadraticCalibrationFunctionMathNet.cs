using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace InternalLogicCalibration
{
    public class QuadraticCalibrationFunctionMathNet : CalibrationFunction
    {
        private Func<double[], double> f;
        private int numFeatures;
        private int numFeaturesExpanded;
        private TransformFunction transform;

        public QuadraticCalibrationFunctionMathNet(TransformFunction transform)
        {
            this.transform = transform;
            numFeatures = transform.numOutputs;
            numFeaturesExpanded = numFeatures + numFeatures * (numFeatures + 1) / 2;
        }

        public override double Predict(double[] t)
        {
            return f(ExpandFeatures(transform.Transform(t)));
        }

        private double[] ExpandFeatures(double[] input)
        {
            double[] outputExpanded = new double[numFeaturesExpanded];
            for (int i = 0; i < numFeatures; i++)
                outputExpanded[i] = input[i];
            int index = numFeatures;
            for (int i = 0; i < numFeatures; i++)
            {
                for (int j = i; j < numFeatures; j++)
                {
                    outputExpanded[index] = input[i] * input[j];
                    index++;
                }
            }
            return outputExpanded;
        }

        public override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            double[][] ok = new double[trainingList.Count()][];
            int k = 0;
            foreach (LabeledDataPoint p in trainingList)
            {
                ok[k] = ExpandFeatures(transform.Transform(p.inputs));
                k++;
            }
            var ok2 = trainingList.Select(b => b.output).ToArray();

            var ye = new Func<double[], double>[numFeaturesExpanded + 1];
            ye[0] = a => 1;
            for (int i = 0; i < numFeaturesExpanded; i++)
            {
                int j = i;
                ye[j + 1] = a => a[j];
            }
            f = Fit.LinearMultiDimFunc(ok, ok2, ye);
            //onOutput("Finished fitting a quadratic"));
        }
    }
}