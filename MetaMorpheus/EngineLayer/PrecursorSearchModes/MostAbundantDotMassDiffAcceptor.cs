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
    /// of shifts, preserving per-shift FDR grouping. Downstream G-PTM-D mod assignment recomputes the
    /// modification from the monoisotopic precursor mass, so it is unaffected by this acceptor.
    /// </remarks>
    public class MostAbundantDotMassDiffAcceptor : MassDiffAcceptor
    {
        private readonly Tolerance Tolerance;
        private readonly AverageResidue Averagine;
        private readonly double[] SortedMassShifts;
        private readonly int[] ShiftNotches;
        public int MaxApexOffsetNeutrons { get; }
        private readonly int[] ApexOffsetsInNeutrons;

        public MostAbundantDotMassDiffAcceptor(string fileNameAddition, IEnumerable<double> acceptableMassShifts,
            Tolerance tol, AverageResidue averagine, int maxApexOffsetNeutrons = 2) : base(fileNameAddition)
        {
            Tolerance = tol;
            Averagine = averagine;
            SortedMassShifts = acceptableMassShifts.OrderBy(Math.Abs).ThenBy(p => p < 0).ToArray();
            ShiftNotches = SortedMassShifts.Select(s => (int)Math.Round(s * NotchScalar)).ToArray();
            MaxApexOffsetNeutrons = maxApexOffsetNeutrons;

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

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            for (int j = 0; j < SortedMassShifts.Length; j++)
            {
                double shiftedMono = peptideMass + SortedMassShifts[j];
                double apex = shiftedMono + ApexOffset(shiftedMono);
                foreach (int k in ApexOffsetsInNeutrons)
                {
                    if (Tolerance.Within(scanPrecursorMass, apex + k * Constants.C13MinusC12))
                    {
                        return ShiftNotches[j];
                    }
                }
            }
            return -1;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
        {
            for (int j = 0; j < SortedMassShifts.Length; j++)
            {
                double shiftedMono = peptideMonoisotopicMass + SortedMassShifts[j];
                double apex = shiftedMono + ApexOffset(shiftedMono);
                foreach (int k in ApexOffsetsInNeutrons)
                {
                    double mass = apex + k * Constants.C13MinusC12;
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
                    double mass = monoApprox - SortedMassShifts[j] - k * Constants.C13MinusC12;
                    yield return new AllowedIntervalWithNotch(Tolerance.GetMinimumValue(mass), Tolerance.GetMaximumValue(mass), ShiftNotches[j]);
                }
            }
        }

        public override string ToString() => FileNameAddition + " mostAbundantDot ±" + MaxApexOffsetNeutrons + " apex " + Tolerance + " " + string.Join(",", SortedMassShifts);

        public override string ToProseString() => Tolerance + " around the most abundant isotopic peak of " + string.Join(", ", SortedMassShifts) + " Da shifts (±" + MaxApexOffsetNeutrons + " isotope apex tolerance)";
    }
}
