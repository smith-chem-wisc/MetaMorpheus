using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace EngineLayer
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

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            double diff = scanPrecursorMass - peptideMass;

            int index = Array.BinarySearch(means, diff);
            if (index >= 0)
                return 0;
            index = ~index;
            // Two options: either it's the index of the first element greater than diff, or len if diff greater than all
            if (index < means.Length && intervals[index].Contains(diff))
                return 0;
            if (index > 0 && intervals[index - 1].Contains(diff))
                return 0;
            return -1;
        }

        public override IEnumerable<Tuple<DoubleRange, int>> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            return intervals.Select(b => new Tuple<DoubleRange, int>(new DoubleRange(peptideMonoisotopicMass + b.Minimum, peptideMonoisotopicMass + b.Maximum), 0));
        }

        #endregion Public Methods

    }
}