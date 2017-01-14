using InternalLogicEngineLayer;
using Spectra;
using System;
using System.Collections.Generic;

namespace InternalLogicTaskLayer
{
	class SingleAbsoluteAroundZeroSearchMode : SearchMode
	{

		#region Private Fields

		readonly double value;

		#endregion Private Fields

		#region Public Constructors

		public SingleAbsoluteAroundZeroSearchMode(string v, double value) : base(v)
		{
			this.value = value;
		}

		#endregion Public Constructors

		#region Public Methods

		public override bool Accepts(double scanPrecursorMass, double peptideMass)
		{
			return Math.Abs(scanPrecursorMass - peptideMass) < value;
		}

		public override IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
		{
			yield return new DoubleRange(-value, value);
		}

		public override string SearchModeString()
		{
			return "SingleAbsoluteAroundZeroSearchMode" + value;
		}

		#endregion Public Methods

	}
}