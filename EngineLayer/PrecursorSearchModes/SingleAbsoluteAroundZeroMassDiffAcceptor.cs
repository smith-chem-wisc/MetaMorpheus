using MzLibUtil;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class SingleAbsoluteAroundZeroSearchMode : MassDiffAcceptor
    {
        private readonly double value;

        public SingleAbsoluteAroundZeroSearchMode(double value) : base(value + "daltonsAroundZero")
        {
            this.value = value;
        }

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            return Math.Abs(scanPrecursorMass - peptideMass) < value ? 0 : -1;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            yield return new AllowedIntervalWithNotch(new DoubleRange(peptideMonoisotopicMass - value, peptideMonoisotopicMass + value), 0);
        }

        public override string ToString()
        {
            return FileNameAddition + " daltonsAroundZero " + value;
        }

        public override string ToProseString()
        {
            return (String.Format("{0:0.000}", value) + " Da around zero");
        }
    }
}