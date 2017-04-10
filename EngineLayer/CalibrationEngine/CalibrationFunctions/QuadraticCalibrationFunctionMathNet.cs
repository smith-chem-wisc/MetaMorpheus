﻿using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Calibration
{
    internal class QuadraticCalibrationFunctionMathNet : CalibrationFunction
    {

        #region Private Fields

        private readonly int numFeatures;
        private readonly int numFeaturesExpanded;
        private Func<double[], double> f;
        private TransformFunction transform;

        #endregion Private Fields

        #region Public Constructors

        public QuadraticCalibrationFunctionMathNet(TransformFunction transform)
        {
            this.transform = transform;
            numFeatures = transform.numOutputs;
            numFeaturesExpanded = numFeatures + numFeatures * (numFeatures + 1) / 2;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("Quadratic");
            sb.Append(" numFeatures: " + numFeatures);
            sb.Append(" numFeaturesExpanded: " + numFeaturesExpanded);
            sb.Append(" transform: " + transform.numOutputs);
            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal override double Predict(double[] t)
        {
            return f(ExpandFeatures(transform.Transform(t)));
        }

        internal override void Train<LabeledDataPoint>(IEnumerable<LabeledDataPoint> trainingList)
        {
            double[][] ok = new double[trainingList.Count()][];
            int k = 0;
            foreach (LabeledDataPoint p in trainingList)
            {
                ok[k] = ExpandFeatures(transform.Transform(p.Inputs));
                k++;
            }
            var ok2 = trainingList.Select(b => b.Label).ToArray();

            var ye = new Func<double[], double>[numFeaturesExpanded + 1];
            ye[0] = a => 1;
            for (int i = 0; i < numFeaturesExpanded; i++)
            {
                int j = i;
                ye[j + 1] = a => a[j];
            }
            f = Fit.LinearMultiDimFunc(ok, ok2, ye);
        }

        #endregion Internal Methods

        #region Private Methods

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

        #endregion Private Methods

    }
}