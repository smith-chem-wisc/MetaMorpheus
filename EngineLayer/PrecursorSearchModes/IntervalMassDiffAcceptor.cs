using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class IntervalMassDiffAcceptor : MassDiffAcceptor
    {
        private readonly List<DoubleRange> Intervals;
        private readonly double[] Means;

        public IntervalMassDiffAcceptor(string fileNameAddition, IEnumerable<DoubleRange> doubleRanges) : base(fileNameAddition)
        {
            Intervals = doubleRanges.OrderBy(b => b.Mean).ToList();
            Means = Intervals.Select(b => b.Mean).ToArray();
        }

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            double diff = scanPrecursorMass - peptideMass;

            int index = Array.BinarySearch(Means, diff);
            if (index >= 0)
                return 0;
            index = ~index;
            // Two options: either it's the index of the first element greater than diff, or len if diff greater than all
            if (index < Means.Length && Intervals[index].Contains(diff))
                return 0;
            if (index > 0 && Intervals[index - 1].Contains(diff))
                return 0;
            return -1;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
        {
            return Intervals.Select(b => new AllowedIntervalWithNotch(new DoubleRange(peptideMonoisotopicMass + b.Minimum, peptideMonoisotopicMass + b.Maximum), 0));
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
        {
            return Intervals.Select(b => new AllowedIntervalWithNotch(new DoubleRange(peptideMonoisotopicMass - b.Maximum, peptideMonoisotopicMass - b.Minimum), 0));
        }

        public override string ToString()
        {
            return FileNameAddition + " interval " + string.Join(",", Intervals);
        }

        public override string ToProseString()
        {
            return ("the mass (Da) interval(s) " + String.Join(", ", Intervals));
        }
    }
}