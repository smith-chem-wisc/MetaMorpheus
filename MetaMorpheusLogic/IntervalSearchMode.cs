using Spectra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MetaMorpheusLogic
{
    public class IntervalSearchMode : SearchMode
    {
        private List<DoubleRange> intervals;

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