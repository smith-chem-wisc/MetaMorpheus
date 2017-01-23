using Spectra;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace InternalLogicEngineLayer
{
    public class IntervalSearchMode : SearchMode
    {

        #region Private Fields

        private readonly List<DoubleRange> intervals;

        #endregion Private Fields

        #region Public Constructors

        public IntervalSearchMode(string fileNameAddition, IEnumerable<DoubleRange> doubleRanges) : base(fileNameAddition)
        {
            intervals = doubleRanges.ToList();
        }

        public IntervalSearchMode(IEnumerable<DoubleRange> doubleRanges) : base("intervals" + string.Join("", doubleRanges.Select(b => "[" + b.Minimum.ToString("F3", CultureInfo.InvariantCulture) + "-" + b.Maximum.ToString("F3", CultureInfo.InvariantCulture) + "]")))
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

        #endregion Public Methods

    }
}