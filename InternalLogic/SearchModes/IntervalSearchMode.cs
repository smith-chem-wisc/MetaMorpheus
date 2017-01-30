using Spectra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace InternalLogicEngineLayer
{
    public class IntervalSearchMode : SearchMode
    {

        #region Private Fields

        private readonly List<DoubleRange> intervals;
        private readonly double[] means;

        #endregion Private Fields

        #region Public Constructors

        public IntervalSearchMode(string fileNameAddition, IEnumerable<DoubleRange> nonOverlappingDoubleRanges) : base(fileNameAddition)
        {
            intervals = nonOverlappingDoubleRanges.OrderBy(b => b.Mean).ToList();
            means = intervals.Select(b => b.Mean).ToArray();
        }

        public IntervalSearchMode(IEnumerable<DoubleRange> doubleRanges) : this("intervals" + string.Join("", doubleRanges.Select(b => "[" + b.Minimum.ToString("F3", CultureInfo.InvariantCulture) + "-" + b.Maximum.ToString("F3", CultureInfo.InvariantCulture) + "]")), doubleRanges)
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override bool Accepts(double scanPrecursorMass, double peptideMass)
        {
            double diff = scanPrecursorMass - peptideMass;

            int index = Array.BinarySearch(means, diff);
            if (index >= 0)
                return true;
            index = ~index;
            // Two options: either it's the index of the first element greater than diff, or len if diff greater than all
            if (index < means.Length)
                if (intervals[index].Contains(diff))
                    return true;
            if (index > 0)
                if (intervals[index - 1].Contains(diff))
                    return true;
            return false;
        }

        public override IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            return intervals.Select(b => new DoubleRange(peptideMonoisotopicMass + b.Minimum, peptideMonoisotopicMass + b.Maximum));
        }

        #endregion Public Methods

    }
}