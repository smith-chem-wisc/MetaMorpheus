using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;

namespace EngineLayer
{
    /// <summary>
    /// Most-abundant-peak analogue of <see cref="DotMassDiffAcceptor"/>, for searches that allow a set
    /// of discrete mass shifts (e.g. G-PTM-D modification masses). For each acceptable shift <c>s</c>,
    /// a candidate is accepted when the observed most-abundant mass matches the theoretical apex of the
    /// shifted candidate—<c>(peptideMono + s) + averagineOffset(peptideMono + s)</c>—within a ±k
    /// neutron apex tolerance (k ∈ [−<see cref="MaxApexOffsetNeutrons"/> .. +]). This composes the
    /// shift set with the most-abundant apex matching and the apex-misprediction notch set described in
    /// <see cref="MostAbundantMassDiffAcceptor"/>.
    /// </summary>
    /// <remarks>
    /// The returned notch identifies the mass <em>shift</em> (same encoding as
    /// <see cref="DotMassDiffAcceptor"/>, round(shiftDa · NotchScalar)); the ±k apex intervals for a
    /// given shift all share that shift's notch, since they represent the same modification differing
    /// only by apex misprediction. <see cref="MassDiffAcceptor.NumNotches"/> therefore equals the number
    /// of shifts, preserving per-shift FDR grouping. Because the notch identifies only the shift, the
    /// apex offset k cannot be read back from it — <see cref="ToObservedMonoisotopicMass"/> re-derives the
    /// (shift, k) pair instead, and downstream G-PTM-D mod assignment MUST go through it: the deconvoluted
    /// monoisotopic mass of an accepted match can be off by whole isotopologues, which would otherwise be
    /// mistaken for a ~1-2 Da modification.
    /// </remarks>
    public class MostAbundantDotMassDiffAcceptor : MassDiffAcceptor
    {
        private readonly Tolerance Tolerance;
        private readonly AverageResidue Averagine;
        private readonly double[] SortedMassShifts;
        private readonly int[] ShiftNotches;
        public int MaxApexOffsetNeutrons { get; }
        private readonly int[] ApexOffsetsInNeutrons;

        /// <summary>
        /// Mass spacing between adjacent isotopologues (the "+1 neutron" step) used to place the apex
        /// notches. Sourced from the deconvolution's
        /// <see cref="MassSpectrometry.DeconvolutionParameters.ExpectedIsotopeSpacing"/> so the acceptor
        /// matches the spacing the envelope was built with — the C-13/C-12 difference for peptides/
        /// proteoforms, but overridable for decoy runs (~0.9444 Da) or non-carbon-dominated polymers.
        /// </summary>
        private readonly double ExpectedIsotopeSpacing;

        public MostAbundantDotMassDiffAcceptor(string fileNameAddition, IEnumerable<double> acceptableMassShifts,
            Tolerance tol, AverageResidue averagine, int maxApexOffsetNeutrons = 2,
            double expectedIsotopeSpacing = Constants.C13MinusC12) : base(fileNameAddition)
        {
            Tolerance = tol;
            Averagine = averagine;
            SortedMassShifts = acceptableMassShifts.OrderBy(Math.Abs).ThenBy(p => p < 0).ToArray();
            ShiftNotches = SortedMassShifts.Select(s => (int)Math.Round(s * NotchScalar)).ToArray();
            MaxApexOffsetNeutrons = maxApexOffsetNeutrons;
            ExpectedIsotopeSpacing = expectedIsotopeSpacing;

            var offsets = new List<int> { 0 };
            for (int k = 1; k <= maxApexOffsetNeutrons; k++) { offsets.Add(-k); offsets.Add(k); }
            ApexOffsetsInNeutrons = offsets.ToArray();

            NumNotches = SortedMassShifts.Length;
        }

        // The averagine most-abundant offset for a monoisotopic mass: the model's diff-to-monoisotopic
        // at the nearest mass bin. (Composes the existing AverageResidue API.) Guards against a large
        // negative shift pushing the shifted mass to/below zero (the averagine model is undefined there),
        // in which case the shift simply won't match rather than indexing the model out of range.
        private double ApexOffset(double monoisotopicMass)
            => monoisotopicMass <= 0 ? 0 : Averagine.GetDiffToMonoisotopic(Averagine.GetMostIntenseMassIndex(monoisotopicMass));

        public override bool MatchesMostAbundantPeak => true;

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            for (int j = 0; j < SortedMassShifts.Length; j++)
            {
                double shiftedMono = peptideMass + SortedMassShifts[j];
                double apex = shiftedMono + ApexOffset(shiftedMono);
                foreach (int k in ApexOffsetsInNeutrons)
                {
                    if (Tolerance.Within(scanPrecursorMass, apex + k * ExpectedIsotopeSpacing))
                    {
                        return ShiftNotches[j];
                    }
                }
            }
            return -1;
        }

        /// <summary>
        /// Undoes the apex offset and the k-neutron apex misprediction for the (shift, k) pair this
        /// acceptor matched on, leaving the observed monoisotopic mass implied by the match — so a
        /// downstream G-PTM-D mod search sees the shift as a chemical mass difference and not as a
        /// mixture of the shift and a whole-isotopologue offset. The apex offset is evaluated at the
        /// SHIFTED monoisotopic mass, exactly as <see cref="Accepts"/> does, because a large shift can
        /// move the candidate into a different averagine apex bin. Falls back to the unshifted apex
        /// model when no (shift, k) pair matches (the caller is then outside this acceptor's decisions).
        /// </summary>
        public override double ToObservedMonoisotopicMass(double observedMostAbundantMass, double theoreticalMonoisotopicMass)
        {
            for (int j = 0; j < SortedMassShifts.Length; j++)
            {
                double shiftedMono = theoreticalMonoisotopicMass + SortedMassShifts[j];
                double apexOffset = ApexOffset(shiftedMono);
                foreach (int k in ApexOffsetsInNeutrons)
                {
                    if (Tolerance.Within(observedMostAbundantMass, shiftedMono + apexOffset + k * ExpectedIsotopeSpacing))
                    {
                        return observedMostAbundantMass - apexOffset - k * ExpectedIsotopeSpacing;
                    }
                }
            }

            double unshiftedApexOffset = ApexOffset(theoreticalMonoisotopicMass);
            double residual = observedMostAbundantMass - (theoreticalMonoisotopicMass + unshiftedApexOffset);
            double nearestK = Math.Round(residual / ExpectedIsotopeSpacing, MidpointRounding.AwayFromZero);
            return observedMostAbundantMass - unshiftedApexOffset - nearestK * ExpectedIsotopeSpacing;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
        {
            for (int j = 0; j < SortedMassShifts.Length; j++)
            {
                double shiftedMono = peptideMonoisotopicMass + SortedMassShifts[j];
                double apex = shiftedMono + ApexOffset(shiftedMono);
                foreach (int k in ApexOffsetsInNeutrons)
                {
                    double mass = apex + k * ExpectedIsotopeSpacing;
                    yield return new AllowedIntervalWithNotch(Tolerance.GetMinimumValue(mass), Tolerance.GetMaximumValue(mass), ShiftNotches[j]);
                }
            }
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double observedMostAbundantMass)
        {
            // Indexed (ModernSearch) path; documented near-boundary caveat as in MostAbundantMassDiffAcceptor.
            // GPTMD uses the theory-driven ClassicSearch path, so this is not its primary path.
            double monoApprox = observedMostAbundantMass - ApexOffset(observedMostAbundantMass);
            for (int j = 0; j < SortedMassShifts.Length; j++)
            {
                foreach (int k in ApexOffsetsInNeutrons)
                {
                    double mass = monoApprox - SortedMassShifts[j] - k * ExpectedIsotopeSpacing;
                    yield return new AllowedIntervalWithNotch(Tolerance.GetMinimumValue(mass), Tolerance.GetMaximumValue(mass), ShiftNotches[j]);
                }
            }
        }

        public override string ToString() => FileNameAddition + " mostAbundantDot ±" + MaxApexOffsetNeutrons + " apex " + Tolerance + " " + string.Join(",", SortedMassShifts);

        public override string ToProseString() => Tolerance + " around the most abundant isotopic peak of " + string.Join(", ", SortedMassShifts) + " Da shifts (±" + MaxApexOffsetNeutrons + " isotope apex tolerance)";
    }
}
