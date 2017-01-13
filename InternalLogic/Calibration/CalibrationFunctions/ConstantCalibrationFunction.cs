using MathNet.Numerics.Statistics;
using System.Collections.Generic;
using System.Linq;

namespace InternalLogicCalibration
{
    public class ConstantCalibrationFunction : CalibrationFunction
    {
        public double a;

        public ConstantCalibrationFunction()
        {
        }

        internal override double Predict(double[] inputs)
        {
            return a;
        }

        internal override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            a = trainingList.Select(b => b.output).Median();
        }
    }
}