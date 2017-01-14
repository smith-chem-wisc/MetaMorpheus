using InternalLogicEngineLayer;
using Spectra;
using System;
using System.Collections.Generic;

namespace InternalLogicTaskLayer
{
    internal class SingleAbsoluteAroundZeroSearchMode : SearchMode
    {

        #region Private Fields

        private string v;
        private double value;

        #endregion Private Fields

        #region Public Constructors

        public SingleAbsoluteAroundZeroSearchMode(string v, double value) : base(v)
        {
            this.v = v;
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
            throw new NotImplementedException();
        }

        public override string SearchModeString()
        {
			return "SingleAbsoluteAroundZeroSearchMode" + value;
        }

        #endregion Public Methods

    }
}