using System.Collections.Generic;

namespace InternalLogicCalibration
{
    public class IdentityCalibrationFunction : CalibrationFunction
    {
        #region Public Constructors

        public IdentityCalibrationFunction()
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override double Predict(double[] inputs)
        {
            return 0;
        }

        public override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
        }

        #endregion Public Methods
    }
}