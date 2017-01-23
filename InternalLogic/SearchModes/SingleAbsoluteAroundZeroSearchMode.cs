using Spectra;
using System;
using System.Collections.Generic;

namespace InternalLogicEngineLayer
{
    public class SingleAbsoluteAroundZeroSearchMode : SearchMode
    {

        #region Private Fields

        private readonly double value;

        #endregion Private Fields

        #region Public Constructors

        public SingleAbsoluteAroundZeroSearchMode(double value) : base(value + "daltonsAroundZero")
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
            yield return new DoubleRange(peptideMonoisotopicMass - value, peptideMonoisotopicMass + value);
        }

        #endregion Public Methods

    }
}