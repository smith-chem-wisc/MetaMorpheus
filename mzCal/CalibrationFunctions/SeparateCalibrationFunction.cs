using System;
using System.Collections.Generic;

namespace mzCal
{
    internal class SeparateCalibrationFunction : CalibrationFunction
    {
        private CalibrationFunction calibrationFunction1;
        private CalibrationFunction calibrationFunction2;

        public SeparateCalibrationFunction(CalibrationFunction calibrationFunction1, CalibrationFunction calibrationFunction2)
        {
            this.calibrationFunction1 = calibrationFunction1;
            this.calibrationFunction2 = calibrationFunction2;
        }

        public override double Predict(double[] inputs)
        {
            if (inputs[0] == 1)
                return calibrationFunction1.Predict(inputs);
            else
                return calibrationFunction2.Predict(inputs);
        }

        public override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            throw new NotImplementedException();
        }
    }
}