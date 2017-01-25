namespace InternalLogicCalibration
{
    internal class SeparateCalibrationFunction : CalibrationFunction
    {

        #region Public Constructors

        public SeparateCalibrationFunction(CalibrationFunction calibrationFunction1, CalibrationFunction calibrationFunction2)
        {
            this.CalibrationFunction1 = calibrationFunction1;
            this.CalibrationFunction2 = calibrationFunction2;
        }

        #endregion Public Constructors

        #region Public Properties

        public CalibrationFunction CalibrationFunction1 { get; private set; }
        public CalibrationFunction CalibrationFunction2 { get; private set; }

        #endregion Public Properties

        #region Internal Methods

        internal override double Predict(double[] t)
        {
            if (t[0] < 0)
                return CalibrationFunction1.Predict(t);
            return CalibrationFunction2.Predict(t);
        }

        #endregion Internal Methods

    }
}