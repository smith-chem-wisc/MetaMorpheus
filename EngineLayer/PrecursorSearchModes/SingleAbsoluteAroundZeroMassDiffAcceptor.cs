using MzLibUtil;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class SingleAbsoluteAroundZeroSearchMode : MassDiffAcceptor
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

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            yield return new AllowedIntervalWithNotch(new DoubleRange(peptideMonoisotopicMass - value, peptideMonoisotopicMass + value), 0);
        }

        public override string ToString()
        {
            return FileNameAddition + " daltonsAroundZero " + value;
        }

        #endregion Public Methods

    }
}