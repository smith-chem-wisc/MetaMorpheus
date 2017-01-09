using System.Collections.Generic;

namespace IndexSearchAndAnalyze
{
    public class IdentityCalibrationFunction : CalibrationFunction
    {
        public IdentityCalibrationFunction()
        {
        }

        public override double Predict(double[] inputs)
        {
            return 0;
        }

        public override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
        }
    }
}