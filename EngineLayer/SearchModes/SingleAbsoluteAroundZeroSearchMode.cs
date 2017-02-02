using Spectra;
using System;
using System.Collections.Generic;

namespace EngineLayer
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

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            return Math.Abs(scanPrecursorMass - peptideMass) < value ? 0 : -1;
        }

        public override IEnumerable<Tuple<DoubleRange, int>> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            yield return new Tuple<DoubleRange, int>(new DoubleRange(peptideMonoisotopicMass - value, peptideMonoisotopicMass + value), 0);
        }

        #endregion Public Methods

    }
}