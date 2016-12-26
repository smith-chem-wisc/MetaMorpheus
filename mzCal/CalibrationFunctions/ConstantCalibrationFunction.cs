using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace mzCal
{
    public class ConstantCalibrationFunction : CalibrationFunction
    {
        public double a;
        private Action<string> onOutput;

        public ConstantCalibrationFunction(Action<string> onOutput)
        {
            this.onOutput = onOutput;
        }

        public override double Predict(double[] inputs)
        {
            return a;
        }

        public override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            a = trainingList.Select(b => b.output).Median();
            onOutput("a = " + a);
        }
    }
}