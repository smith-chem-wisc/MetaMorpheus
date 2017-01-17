using Spectra;
using System;
using System.Collections.Generic;

namespace InternalLogicEngineLayer
{
    public class SinglePpmAroundZeroSearchMode : SearchMode
    {

        #region Private Fields

        private readonly double ppmTolerance;

        #endregion Private Fields

        #region Public Constructors

        public SinglePpmAroundZeroSearchMode(string v1, double ppmTolerance) : base(v1)
        {
            this.ppmTolerance = ppmTolerance;
        }

        #endregion Public Constructors

        #region Public Methods

        public override bool Accepts(double scanPrecursorMass, double peptideMass)
        {
            return Math.Abs((scanPrecursorMass - peptideMass) / (peptideMass) * 1e6) < ppmTolerance;
        }

        public override IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            var diff = ppmTolerance / 1e6 * peptideMonoisotopicMass;
            yield return new DoubleRange(peptideMonoisotopicMass - diff, peptideMonoisotopicMass + diff);
        }

        public override string SearchModeString()
        {
            return "SinglePpmAroundZeroSearchMode" + ppmTolerance;
        }

        #endregion Public Methods

    }
}