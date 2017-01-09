using MathNet.Numerics.Statistics;
using System.Collections.Generic;
using System.Linq;

namespace IndexSearchAndAnalyze
{
    public class ConstantCalibrationFunction : CalibrationFunction
    {
        public double a;

        public ConstantCalibrationFunction()
        {
        }

        public override double Predict(double[] inputs)
        {
            return a;
        }

        public override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            a = trainingList.Select(b => b.output).Median();
        }
    }
}