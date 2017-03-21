using System;
using System.Collections.Generic;

namespace EngineLayer.Calibration
{
    public abstract class CalibrationFunction
    {

        #region Internal Methods

        internal abstract double Predict(double[] t);

        internal abstract void Train<LabeledDataPoint>(IEnumerable<LabeledDataPoint> trainingList)
            where LabeledDataPoint : IHasInputsAndOutputs;

        internal double GetMSE<LabeledDataPoint>(IEnumerable<LabeledDataPoint> pointList)
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

        #endregion Internal Methods

    }
}