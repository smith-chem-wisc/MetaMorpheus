using System;
using System.Collections.Generic;

namespace mzCal
{
    public class IdentityCalibrationFunction : CalibrationFunction
    {
        private Action<string> onOutput;

        public IdentityCalibrationFunction(Action<string> onOutput)
        {
            this.onOutput = onOutput;
        }

        public override double Predict(double[] inputs)
        {
            return 0;
        }

        public override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            onOutput("Sucessfully trained IdentityCalibrationFunction");
        }
    }
}