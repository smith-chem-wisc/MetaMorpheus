using Spectra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace InternalLogicEngineLayer
{
    public class IntervalSearchMode : SearchMode
    {
        #region Private Fields

        private List<DoubleRange> intervals;

        #endregion Private Fields

        #region Public Constructors

        public IntervalSearchMode(string fileNameAddition, IEnumerable<DoubleRange> doubleRanges) : base(fileNameAddition)
        {
            intervals = doubleRanges.ToList();
        }

        #endregion Public Constructors

        #region Public Methods

        public override bool Accepts(double scanPrecursorMass, double peptideMass)
        {
            throw new NotImplementedException();
        }

        #endregion Public Methods

        #region Internal Methods

        internal override IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            return intervals.Select(b => new DoubleRange(peptideMonoisotopicMass + b.Minimum, peptideMonoisotopicMass + b.Maximum));
        }

        internal override string SearchModeString()
        {
            return "Intervals allowed: " + string.Join(",", intervals);
        }

        #endregion Internal Methods
    }
}