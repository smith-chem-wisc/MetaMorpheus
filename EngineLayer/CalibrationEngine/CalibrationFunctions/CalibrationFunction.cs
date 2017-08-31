using System;
using System.Collections.Generic;

namespace EngineLayer.Calibration
{
    public abstract class CalibrationFunction
    {
        #region Public Methods

        public abstract double Predict(double[] t);

        public abstract void Train<LabeledDataPoint>(IEnumerable<LabeledDataPoint> trainingList)
            where LabeledDataPoint : IHasInputsAndOutputs;

        public double GetMSE<LabeledDataPoint>(IEnumerable<LabeledDataPoint> pointList)
            where LabeledDataPoint : IHasInputsAndOutputs
        {
            double mse = 0;
            int count = 0;
            foreach (LabeledDataPoint p in pointList)
            {
                mse += Math.Pow(Predict(p.Inputs) - p.Label, 2);
                count++;
            }
            return count == 0 ? 0 : mse / count;
        }

        #endregion Public Methods
    }
}