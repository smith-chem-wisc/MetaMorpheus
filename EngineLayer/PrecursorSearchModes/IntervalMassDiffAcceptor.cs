﻿using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class IntervalMassDiffAcceptor : MassDiffAcceptor
    {
        private readonly List<DoubleRange> intervals;
        private readonly double[] means;

        public IntervalMassDiffAcceptor(string fileNameAddition, IEnumerable<DoubleRange> doubleRanges) : base(fileNameAddition)
        {
            intervals = doubleRanges.OrderBy(b => b.Mean).ToList();
            means = intervals.Select(b => b.Mean).ToArray();
        }

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

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            return intervals.Select(b => new AllowedIntervalWithNotch(new DoubleRange(peptideMonoisotopicMass + b.Minimum, peptideMonoisotopicMass + b.Maximum), 0));
        }

        public override string ToString()
        {
            return FileNameAddition + " interval " + string.Join(",", intervals);
        }

        public override string ToProseString()
        {
            return ("the mass (Da) interval(s) " + String.Join(", ", intervals));
        }
    }
}