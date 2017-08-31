using System.Collections.Generic;
using System.Text;

namespace EngineLayer.Calibration
{
    public class IdentityCalibrationFunction : CalibrationFunction
    {
        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("Identity");
            return sb.ToString();
        }

        public override double Predict(double[] t)
        {
            return 0;
        }

        public override void Train<LabeledDataPoint>(IEnumerable<LabeledDataPoint> trainingList)
        {
            // Does nothing because this is the identity calibration function
        }

        #endregion Public Methods
    }
}