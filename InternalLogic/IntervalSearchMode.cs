using Spectra;
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
            foreach (var huh in intervals)
            {
                if (huh.Contains(scanPrecursorMass - peptideMass))
                    return true;
            }
            return false;
        }

        public override IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            return intervals.Select(b => new DoubleRange(peptideMonoisotopicMass + b.Minimum, peptideMonoisotopicMass + b.Maximum));
        }

        public override string SearchModeString()
        {
            return "Intervals allowed: " + string.Join(",", intervals);
        }

        #endregion Public Methods

    }
}