using MathNet.Numerics.Statistics;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Calibration
{
    public class ConstantCalibrationFunction : CalibrationFunction
    {

        #region Public Fields

        public double a;

        #endregion Public Fields

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("Constant");
            sb.Append(" a = " + a);
            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal override double Predict(double[] t)
        {
            return a;
        }

        internal override void Train<LabeledDataPoint>(IEnumerable<LabeledDataPoint> trainingList)
        {
            a = trainingList.Select(b => b.Label).Median();
        }

        #endregion Internal Methods

    }
}