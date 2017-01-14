namespace InternalLogicCalibration
{
	class SeparateCalibrationFunction : CalibrationFunction
	{
		readonly CalibrationFunction calibrationFunction1;
		readonly CalibrationFunction calibrationFunction2;

		public SeparateCalibrationFunction(CalibrationFunction calibrationFunction1, CalibrationFunction calibrationFunction2)
		{
			this.calibrationFunction1 = calibrationFunction1;
			this.calibrationFunction2 = calibrationFunction2;
		}

		internal override double Predict(double[] t)
		{
			if (t[0] < 0)
				return calibrationFunction1.Predict(t);
			return calibrationFunction2.Predict(t);
		}
	}
}