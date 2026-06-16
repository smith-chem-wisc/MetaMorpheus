using System.Collections.Generic;
using MassSpectrometry;
using MzLibUtil;

namespace EngineLayer
{
    /// <summary>
    /// Selects theoretical candidates by matching the observed most-abundant isotopic peak mass
    /// against each candidate's theoretical most-abundant mass (Strategy B). The theoretical
    /// most-abundant mass is the candidate's exact monoisotopic mass plus an averagine-derived,
    /// mass-dependent offset (roughly an integer number of neutrons). Because the same physical
    /// (most-abundant) peak is compared on both the observed and theoretical sides, this avoids the
    /// monoisotopic off-by-N errors that arise when the true monoisotopic peak is undetectable.
    ///
    /// Used when <see cref="CommonParameters.PrecursorMassMatchMode"/> is
    /// <see cref="PrecursorMassMatchMode.MostAbundant"/>. A single notch (0) is returned: this
    /// acceptor replaces the missed-monoisotopic notches, which exist precisely to paper over the
    /// off-by-N problem this mode solves directly.
    /// </summary>
    /// <remarks>
    /// The scan-side mass supplied to <see cref="Accepts"/> /
    /// <see cref="GetAllowedPrecursorMassIntervalsFromObservedMass"/> is expected to be the envelope's
    /// most-abundant observed neutral mass (set on the scan as PrecursorMassToMatch). Isotopically
    /// unresolved species would instead supply an average mass and require
    /// <see cref="AverageResidue.GetAverageOffset"/>; that path is deferred to v2 along with the
    /// unresolved charge-determination algorithm, so all envelopes are currently resolved.
    /// </remarks>
    public class MostAbundantMassDiffAcceptor : MassDiffAcceptor
    {
        private readonly Tolerance Tolerance;
        private readonly AverageResidue Averagine;

        public MostAbundantMassDiffAcceptor(string fileNameAddition, Tolerance tol, AverageResidue averagine)
            : base(fileNameAddition)
        {
            Tolerance = tol;
            Averagine = averagine;
            NumNotches = 1;
        }

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            double theoreticalMostAbundantMass = peptideMass + Averagine.GetMostAbundantOffset(peptideMass);
            return Tolerance.Within(scanPrecursorMass, theoreticalMostAbundantMass) ? 0 : -1;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
        {
            double mass = peptideMonoisotopicMass + Averagine.GetMostAbundantOffset(peptideMonoisotopicMass);
            yield return new AllowedIntervalWithNotch(Tolerance.GetMinimumValue(mass), Tolerance.GetMaximumValue(mass), 0);
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
        {
            // Indexed (ModernSearch) path: convert the observed most-abundant mass back to a candidate
            // monoisotopic window by subtracting the offset evaluated near the observed mass. Top-down
            // uses the theory-driven ClassicSearch path (GetAllowed...FromTheoreticalMass), which is
            // exact; this observed-side conversion carries a small near-boundary ambiguity and is only
            // exercised by indexed bottom-up search, which this mode does not target.
            double mass = peptideMonoisotopicMass - Averagine.GetMostAbundantOffset(peptideMonoisotopicMass);
            yield return new AllowedIntervalWithNotch(Tolerance.GetMinimumValue(mass), Tolerance.GetMaximumValue(mass), 0);
        }

        public override string ToString()
        {
            return FileNameAddition + " mostAbundant " + Tolerance;
        }

        public override string ToProseString()
        {
            return Tolerance + " around the most abundant isotopic peak";
        }
    }
}
