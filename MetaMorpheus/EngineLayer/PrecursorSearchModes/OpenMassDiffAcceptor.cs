using MzLibUtil;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class OpenSearchMode : MassDiffAcceptor
    {
        public OpenSearchMode() : base("OpenSearch")
        {
        }

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            return 0;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
        {
            yield return new AllowedIntervalWithNotch(new DoubleRange(Double.NegativeInfinity, Double.PositiveInfinity), 0);
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
        {
            yield return new AllowedIntervalWithNotch(new DoubleRange(Double.NegativeInfinity, Double.PositiveInfinity), 0);
        }

        public override string ToProseString()
        {
            return ("unbounded");
        }

        public override string ToString()
        {
            return FileNameAddition + " OpenSearch";
        }
    }
}