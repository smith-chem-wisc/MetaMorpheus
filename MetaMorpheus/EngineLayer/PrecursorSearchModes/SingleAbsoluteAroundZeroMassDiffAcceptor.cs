using MzLibUtil;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class SingleAbsoluteAroundZeroSearchMode : MassDiffAcceptor
    {
        private readonly double Value;

        public SingleAbsoluteAroundZeroSearchMode(double value) : base(value + "daltonsAroundZero")
        {
            this.Value = value;
        }

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            return Math.Abs(scanPrecursorMass - peptideMass) < Value ? 0 : -1;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
        {
            yield return new AllowedIntervalWithNotch(new DoubleRange(peptideMonoisotopicMass - Value, peptideMonoisotopicMass + Value), 0);
        }
        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
        {
            yield return new AllowedIntervalWithNotch(new DoubleRange(peptideMonoisotopicMass - Value, peptideMonoisotopicMass + Value), 0);
        }

        public override string ToString()
        {
            return FileNameAddition + " daltonsAroundZero " + Value;
        }

        public override string ToProseString()
        {
            return (String.Format("{0:0.000}", Value) + " Da around zero");
        }
    }
}