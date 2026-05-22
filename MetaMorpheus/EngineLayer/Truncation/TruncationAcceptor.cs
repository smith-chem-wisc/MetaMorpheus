using System;
using System.Collections.Generic;
using MzLibUtil;

namespace EngineLayer.Truncation
{
    /// <summary>
    /// <see cref="MassDiffAcceptor"/> for the truncation search Pass 2 precursor filter (decision #7).
    /// A parent of theoretical mass M_theo is a candidate for an observed precursor mass M_obs iff
    /// <c>0 &lt; M_obs &lt;= M_theo + tolerance</c> — the parent must be at least as heavy as the scan,
    /// allowing M_obs to exceed M_theo by no more than the precursor tolerance. Single category, no notches.
    ///
    /// Phase 0 stub: the acceptance logic is implemented in Phase 1.
    /// </summary>
    public class TruncationAcceptor : MassDiffAcceptor
    {
        private readonly Tolerance _precursorMassTolerance;

        public TruncationAcceptor(Tolerance precursorMassTolerance) : base("truncation")
        {
            _precursorMassTolerance = precursorMassTolerance;
        }

        // TODO Phase 1: implement 0 < M_obs <= M_theo + tolerance (single notch 0).
        public override int Accepts(double scanPrecursorMass, double peptideMass) =>
            throw new NotImplementedException();

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass) =>
            throw new NotImplementedException();

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass) =>
            throw new NotImplementedException();

        public override string ToProseString() =>
            $"truncation precursor acceptor (0 < observed <= theoretical + {_precursorMassTolerance})";
    }
}
