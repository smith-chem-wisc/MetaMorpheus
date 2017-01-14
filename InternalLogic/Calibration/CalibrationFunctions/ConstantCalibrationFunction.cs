using MathNet.Numerics.Statistics;
using System.Collections.Generic;
using System.Linq;

namespace InternalLogicCalibration
{
	public class ConstantCalibrationFunction : CalibrationFunction
	{
		public double a;

		internal override double Predict(double[] t)
		{
			return a;
		}

		internal void Train(IEnumerable<LabeledDataPoint> trainingList)
		{
			a = trainingList.Select(b => b.output).Median();
		}
	}
}