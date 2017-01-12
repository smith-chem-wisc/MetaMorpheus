using System;
using System.Collections.Generic;

namespace InternalLogicCalibration
{
    internal class ByHandCalibrationFunction : CalibrationFunction
    {
        #region Private Fields

        private Action<string> onOutput;

        #endregion Private Fields

        #region Public Constructors

        public ByHandCalibrationFunction(Action<string> onOutput)
        {
            this.onOutput = onOutput;
        }

        #endregion Public Constructors

        #region Public Methods

        internal override double Predict(double[] t)
        {
            return -t[1] / 200000;
        }

        internal override void Train(IEnumerable<LabeledDataPoint> trainingList)
        {
            onOutput("Sucessfully trained ByHandCalibrationFunction");
        }

        #endregion Public Methods
    }
}