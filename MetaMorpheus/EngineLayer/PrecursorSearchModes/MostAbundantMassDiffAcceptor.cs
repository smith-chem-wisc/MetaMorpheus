using System;
using System.Collections.Generic;
using Chemistry;
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
    /// <see cref="PrecursorMassMatchMode.MostAbundant"/>.
    /// </summary>
    /// <remarks>
    /// The averagine offset predicts which isotopologue is tallest, but for real proteoforms the
    /// experimental apex can differ from the averagine prediction by ±1–2 neutrons (validated on
    /// yeast top-down data: a single tight-ppm point at the predicted apex under-identifies badly,
    /// because one neutron is ~67 ppm at 15 kDa). To tolerate that apex misprediction while keeping
    /// each match at tight ppm (and FDR controlled per-notch), this acceptor emits a small set of
    /// notches at <c>apex + k·C13</c> for k in [−<see cref="MaxApexOffsetNeutrons"/> ..
    /// +<see cref="MaxApexOffsetNeutrons"/>]. k = 0 is the confident on-apex match; nonzero k are the
    /// apex-misprediction cases. Notch integers follow the existing DotMassDiffAcceptor convention
    /// (round(shiftDa · NotchScalar)) so they never collide with the −1 "not accepted" sentinel.
    ///
    /// The scan-side mass is the envelope's most-abundant observed neutral mass (PrecursorMassToMatch).
    /// Isotopically unresolved species would supply an average mass and need
    /// <see cref="AverageResidue.GetAverageOffset"/>; that path is deferred to v2 with the unresolved
    /// charge-determination algorithm, so all envelopes are currently resolved.
    /// </remarks>
    public class MostAbundantMassDiffAcceptor : MassDiffAcceptor
    {
        private readonly Tolerance Tolerance;
        private readonly AverageResidue Averagine;

        /// <summary>Maximum apex misprediction, in neutrons, tolerated on either side of the averagine-predicted apex.</summary>
        public int MaxApexOffsetNeutrons { get; }

        // k values ordered by ascending |k| (then sign) so the on-apex (k = 0) match is preferred.
        private readonly int[] ApexOffsetsInNeutrons;

        public MostAbundantMassDiffAcceptor(string fileNameAddition, Tolerance tol, AverageResidue averagine, int maxApexOffsetNeutrons = 2)
            : base(fileNameAddition)
        {
            Tolerance = tol;
            Averagine = averagine;
            MaxApexOffsetNeutrons = maxApexOffsetNeutrons;

            var offsets = new List<int> { 0 };
            for (int k = 1; k <= maxApexOffsetNeutrons; k++)
            {
                offsets.Add(-k);
                offsets.Add(k);
            }
            ApexOffsetsInNeutrons = offsets.ToArray();
            NumNotches = ApexOffsetsInNeutrons.Length;
        }

        private static int NotchFor(int k) => (int)Math.Round(k * Constants.C13MinusC12 * NotchScalar);

        // The averagine most-abundant offset for a monoisotopic mass: the model's diff-to-monoisotopic
        // at the nearest mass bin. (Composes the existing AverageResidue API rather than relying on a
        // dedicated mzLib method.)
        private double ApexOffset(double monoisotopicMass)
            => Averagine.GetDiffToMonoisotopic(Averagine.GetMostIntenseMassIndex(monoisotopicMass));

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            double apex = peptideMass + ApexOffset(peptideMass);
            foreach (int k in ApexOffsetsInNeutrons) // ordered to prefer k = 0
            {
                if (Tolerance.Within(scanPrecursorMass, apex + k * Constants.C13MinusC12))
                {
                    return NotchFor(k);
                }
            }
            return -1;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
        {
            double apex = peptideMonoisotopicMass + ApexOffset(peptideMonoisotopicMass);
            foreach (int k in ApexOffsetsInNeutrons)
            {
                double mass = apex + k * Constants.C13MinusC12;
                yield return new AllowedIntervalWithNotch(Tolerance.GetMinimumValue(mass), Tolerance.GetMaximumValue(mass), NotchFor(k));
            }
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
        {
            // Indexed (ModernSearch) path: convert the observed most-abundant mass back to candidate
            // monoisotopic windows by subtracting the offset evaluated near the observed mass, for each
            // apex-offset notch. Top-down uses the exact theory-driven ClassicSearch path; this
            // observed-side conversion carries a small near-boundary ambiguity and is only exercised by
            // indexed bottom-up search, which this mode does not target.
            double monoApprox = peptideMonoisotopicMass - ApexOffset(peptideMonoisotopicMass);
            foreach (int k in ApexOffsetsInNeutrons)
            {
                double mass = monoApprox - k * Constants.C13MinusC12;
                yield return new AllowedIntervalWithNotch(Tolerance.GetMinimumValue(mass), Tolerance.GetMaximumValue(mass), NotchFor(k));
            }
        }

        public override string ToString()
        {
            return FileNameAddition + " mostAbundant ±" + MaxApexOffsetNeutrons + " apex " + Tolerance;
        }

        public override string ToProseString()
        {
            return Tolerance + " around the most abundant isotopic peak (±" + MaxApexOffsetNeutrons + " isotope apex tolerance)";
        }
    }
}
