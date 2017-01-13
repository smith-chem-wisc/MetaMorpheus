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

        internal override double Predict(double[] inputs)
        {
            return 0;
        }

        internal override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
        }

        #endregion Public Methods
    }
}