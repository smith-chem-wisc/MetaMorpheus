using System;
using System.Collections.Generic;
using MzLibUtil;

namespace EngineLayer.Truncation
{
    /// <summary>
    /// <see cref="MassDiffAcceptor"/> for the truncation search Pass 2 precursor filter (decision #7).
    /// A parent of theoretical mass <c>M_theo</c> is a candidate for an observed precursor mass
    /// <c>M_obs</c> iff <c>0 &lt; M_obs &lt;= M_theo + tolerance</c> — the parent must be at least as
    /// heavy as the scan, allowing M_obs to exceed M_theo by no more than the precursor tolerance
    /// (which absorbs measurement error and the intact-equality case). Single category, single notch (0).
    /// The tolerance band is evaluated on the theoretical mass, consistent with the rest of MetaMorpheus.
    /// </summary>
    public class TruncationAcceptor : MassDiffAcceptor
    {
        private readonly Tolerance _precursorMassTolerance;

        public TruncationAcceptor(Tolerance precursorMassTolerance) : base("truncation")
        {
            _precursorMassTolerance = precursorMassTolerance;
            // base sets NumNotches = 1; this acceptor has a single category, notch 0.
        }

        /// <summary>Returns 0 (accepted, notch 0) iff <c>0 &lt; scanPrecursorMass &lt;= peptideMass + tolerance</c>, else -1.</summary>
        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            if (scanPrecursorMass <= 0)
            {
                return -1;
            }

            // peptideMass == M_theo (the parent). Accept when the observed precursor does not exceed
            // the parent mass by more than the tolerance (i.e. the parent is heavy enough to be truncated
            // down to the observed mass).
            return scanPrecursorMass <= _precursorMassTolerance.GetMaximumValue(peptideMass) ? 0 : -1;
        }

        /// <summary>
        /// Given an observed precursor mass, the acceptable theoretical (parent) masses are
        /// <c>[M_obs - tolerance, +inf)</c> — any parent at least as heavy as the scan (within tolerance).
        /// Used for binary-search pruning of fragment-index bins; <see cref="Accepts"/> enforces the exact rule.
        /// </summary>
        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
        {
            yield return new AllowedIntervalWithNotch(
                _precursorMassTolerance.GetMinimumValue(peptideMonoisotopicMass),
                double.PositiveInfinity,
                0);
        }

        /// <summary>
        /// Given a theoretical (parent) mass, the acceptable observed masses are
        /// <c>(0, M_theo + tolerance]</c>.
        /// </summary>
        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
        {
            yield return new AllowedIntervalWithNotch(
                0,
                _precursorMassTolerance.GetMaximumValue(peptideMonoisotopicMass),
                0);
        }

        public override string ToProseString() =>
            $"truncation precursor acceptor (0 < observed <= theoretical + {_precursorMassTolerance})";
    }
}
