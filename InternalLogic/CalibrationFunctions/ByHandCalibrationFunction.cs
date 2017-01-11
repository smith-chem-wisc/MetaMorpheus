using System;
using System.Collections.Generic;

namespace InternalLogic
{
    internal class ByHandCalibrationFunction : CalibrationFunction
    {
        private Action<string> onOutput;

        public ByHandCalibrationFunction(Action<string> onOutput)
        {
            this.onOutput = onOutput;
        }

        public override double Predict(double[] t)
        {
            return -t[1] / 200000;
        }

        public override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            onOutput("Sucessfully trained ByHandCalibrationFunction");
        }
    }
}