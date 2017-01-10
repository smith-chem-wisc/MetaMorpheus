using Spectra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace IndexSearchAndAnalyze
{
    public class IntervalSearchMode : SearchMode
    {
        List<DoubleRange> intervals;
        public IntervalSearchMode(string fileNameAddition, IEnumerable<DoubleRange> doubleRanges) : base(fileNameAddition)
        {
            intervals = doubleRanges.ToList();
        }

        public override bool Accepts(double scanPrecursorMass, double peptideMass)
        {
            throw new NotImplementedException();
        }
    }
}