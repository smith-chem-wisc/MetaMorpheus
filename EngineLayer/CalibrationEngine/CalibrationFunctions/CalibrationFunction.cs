using System;
using System.Collections.Generic;

namespace EngineLayer.Calibration
{
    public abstract class CalibrationFunction
    {

        #region Internal Methods

        internal abstract double Predict(double[] t);

        internal double getMSE(IEnumerable<LabeledDataPoint> pointList)
        {
            double mse = 0;
            int count = 0;
            foreach (LabeledDataPoint p in pointList)
            {
                mse += Math.Pow(Predict(p.inputs) - p.output, 2);
                count++;
            }
            return count == 0 ? 0 : mse / count;
        }

        #endregion Internal Methods

    }
}