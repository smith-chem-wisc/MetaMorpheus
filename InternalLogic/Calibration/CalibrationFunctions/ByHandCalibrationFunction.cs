using System;

namespace InternalLogicCalibration
{
    internal class ByHandCalibrationFunction : CalibrationFunction
    {
        #region Private Fields

        private readonly Action<string> onOutput;

        #endregion Private Fields

        #region Public Constructors

        public ByHandCalibrationFunction(Action<string> onOutput)
        {
            this.onOutput = onOutput;
        }

        #endregion Public Constructors

        #region Internal Methods

        internal override double Predict(double[] t)
        {
            return -t[1] / 200000;
        }

        internal void Train()
        {
            onOutput("Sucessfully trained ByHandCalibrationFunction");
        }

        #endregion Internal Methods
    }
}