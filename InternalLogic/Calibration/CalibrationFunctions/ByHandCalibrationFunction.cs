using System;

namespace InternalLogicCalibration
{
	class ByHandCalibrationFunction : CalibrationFunction
	{
		#region Private Fields

		readonly Action<string> onOutput;

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

		internal void Train()
		{
			onOutput("Sucessfully trained ByHandCalibrationFunction");
		}

		#endregion Public Methods
	}
}