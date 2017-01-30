using Spectra;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class SinglePpmAroundZeroSearchMode : SearchMode
    {

        #region Private Fields

        private readonly double ppmTolerance;

        #endregion Private Fields

        #region Public Constructors

        public SinglePpmAroundZeroSearchMode(double ppmTolerance) : base(ppmTolerance + "ppmAroundZero")
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

        #endregion Public Methods

    }
}