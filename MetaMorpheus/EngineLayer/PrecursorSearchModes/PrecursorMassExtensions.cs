using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    /// <summary>
    /// Reads the precursor mass of a scan or a spectral match the way the CURRENT SEARCH means it.
    /// <para>
    /// <see cref="Ms2ScanWithSpecificMass"/> and <see cref="SpectralMatch"/> store only OBSERVATIONS —
    /// the deconvoluted monoisotopic mass and the most-abundant (tallest) isotopologue mass. Neither
    /// property changes meaning with the search being run, so reading
    /// <see cref="Ms2ScanWithSpecificMass.PrecursorMass"/> always gives the monoisotopic mass and can
    /// never silently return something else. The search type lives in
    /// <see cref="CommonParameters.PrecursorMassMatchMode"/>, which every engine already has, so these
    /// extensions take the parameters and do the interpreting in one place.
    /// </para>
    /// <para>
    /// Two different questions are answered here, and using the wrong one is a bug:
    /// <list type="bullet">
    /// <item><see cref="GetPrecursorMassForSearch(Ms2ScanWithSpecificMass, CommonParameters)"/> — which
    /// observed mass do we SELECT candidates with.</item>
    /// <item><see cref="GetObservedMonoisotopicMass"/> / <see cref="MatchesCandidateMass"/> — how do we read
    /// a mass DIFFERENCE against a candidate without mistaking an isotope-assignment offset for chemistry.</item>
    /// </list>
    /// </para>
    /// </summary>
    public static class PrecursorMassExtensions
    {
        // Allocating an Averagine builds its mass table, so keep one for the (rare) case where a search
        // runs without precursor deconvolution parameters.
        private static readonly AverageResidue DefaultAveragine = new Averagine();

        /// <summary>
        /// The observed precursor mass used to select theoretical candidates: the monoisotopic mass in the
        /// default mode, the most-abundant isotopologue mass in
        /// <see cref="PrecursorMassMatchMode.MostAbundant"/>. Falls back to the monoisotopic mass when no
        /// most-abundant peak was observed — no deconvoluted envelope (a scan-header precursor, or a
        /// neutral mass read from a pre-deconvoluted file). A most-abundant search cannot match on a peak
        /// that was never seen.
        /// </summary>
        public static double GetPrecursorMassForSearch(this Ms2ScanWithSpecificMass scan, CommonParameters commonParameters)
            => UsesMostAbundantPeak(commonParameters) && scan.PrecursorMostAbundantMass > 0
                ? scan.PrecursorMostAbundantMass
                : scan.PrecursorMass;

        /// <inheritdoc cref="GetPrecursorMassForSearch(Ms2ScanWithSpecificMass, CommonParameters)"/>
        public static double GetPrecursorMassForSearch(this SpectralMatch psm, CommonParameters commonParameters)
            => UsesMostAbundantPeak(commonParameters) && psm.ScanPrecursorMostAbundantMass > 0
                ? psm.ScanPrecursorMostAbundantMass
                : psm.ScanPrecursorMass;

        /// <summary>
        /// The observed MONOISOTOPIC mass implied by this match, with the isotope-level offsets a
        /// most-abundant search allowed removed. Subtracting the candidate's theoretical monoisotopic mass
        /// from the result gives the CHEMICAL mass difference — what localization localizes, what a
        /// histogram bins, and what calibration reads as instrument error.
        /// <para>
        /// A most-abundant search accepts candidates whose deconvoluted monoisotopic peak is off by whole
        /// isotopologues (the apex notch set) — those PSMs are precisely the ones the mode exists to
        /// rescue. For them, (ScanPrecursorMass - theoretical) carries a k-neutron offset that is an
        /// isotope-assignment artifact, not chemistry. Read raw, it becomes a phantom ~1-2 Da modification.
        /// </para>
        /// <para>
        /// The k the search matched on is recovered by rounding the residual beyond the predicted apex to
        /// the nearest whole isotopologue. That is exact BECAUSE most-abundant mode always searches with the
        /// zero-centred apex acceptor (<see cref="SearchTask"/>.GetMassDiffAcceptor overrides the
        /// MassDiffAcceptorType), so the residual is k neutrons plus a within-tolerance ppm error and
        /// nothing else. It is NOT valid for a search carrying discrete mass shifts, where an unknown
        /// modification mass would be folded into k — G-PTM-D must use
        /// <see cref="MatchesCandidateMass"/> instead, which tests each candidate mass on its own.
        /// </para>
        /// Returns <see cref="SpectralMatch.ScanPrecursorMass"/> unchanged in the default monoisotopic mode,
        /// so baseline behaviour is untouched.
        /// </summary>
        public static double GetObservedMonoisotopicMass(this SpectralMatch psm, double theoreticalMonoisotopicMass, CommonParameters commonParameters)
        {
            if (!UsesMostAbundantPeak(commonParameters) || psm.ScanPrecursorMostAbundantMass <= 0 || theoreticalMonoisotopicMass <= 0)
            {
                return psm.ScanPrecursorMass;
            }

            double observedApex = psm.ScanPrecursorMostAbundantMass;
            double apexOffset = ApexOffset(commonParameters, theoreticalMonoisotopicMass);
            double isotopeSpacing = IsotopeSpacing(commonParameters);

            double residual = observedApex - (theoreticalMonoisotopicMass + apexOffset);
            // AwayFromZero so an exact half-isotopologue residual is parity-independent (the default
            // banker's rounding would snap to the even k).
            double neutronOffsets = Math.Round(residual / isotopeSpacing, MidpointRounding.AwayFromZero);

            return observedApex - apexOffset - neutronOffsets * isotopeSpacing;
        }

        /// <summary>
        /// Whether the observed precursor supports a candidate of the given monoisotopic mass, judged the
        /// way the current search judges it: monoisotopic-to-monoisotopic in the default mode, or
        /// apex-to-apex within a ±k isotopologue apex tolerance in most-abundant mode.
        /// <para>
        /// This is the form a search carrying discrete mass shifts needs (G-PTM-D): each candidate mass
        /// (peptide + modification) is tested on its own, so the modification mass can never be confused
        /// with the k-neutron apex misprediction. It mirrors the acceptance test in
        /// <see cref="MostAbundantDotMassDiffAcceptor.Accepts"/> — both call the same apex model — but needs
        /// only the deconvolution parameters, not the acceptor's shift set.
        /// </para>
        /// </summary>
        public static bool MatchesCandidateMass(this SpectralMatch psm, double candidateMonoisotopicMass, Tolerance tolerance, CommonParameters commonParameters)
        {
            if (!UsesMostAbundantPeak(commonParameters) || psm.ScanPrecursorMostAbundantMass <= 0 || candidateMonoisotopicMass <= 0)
            {
                return tolerance.Within(psm.ScanPrecursorMass, candidateMonoisotopicMass);
            }

            double apex = candidateMonoisotopicMass + ApexOffset(commonParameters, candidateMonoisotopicMass);
            double isotopeSpacing = IsotopeSpacing(commonParameters);

            for (int k = -MostAbundantMassDiffAcceptor.DefaultMaxApexOffsetNeutrons; k <= MostAbundantMassDiffAcceptor.DefaultMaxApexOffsetNeutrons; k++)
            {
                if (tolerance.Within(psm.ScanPrecursorMostAbundantMass, apex + k * isotopeSpacing))
                {
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// A warning for a most-abundant search that will silently match some scans on their monoisotopic
        /// mass, because no most-abundant peak was observed for them. Returns <c>null</c> when there is
        /// nothing to report: a monoisotopic search (where the fallback does not exist), an empty scan set,
        /// or a scan set in which every precursor has an observed apex.
        /// <para>
        /// The fallback in <see cref="GetPrecursorMassForSearch(Ms2ScanWithSpecificMass, CommonParameters)"/>
        /// is deliberate — a most-abundant search cannot match on a peak that was never seen — but it is
        /// invisible at the call site, and a run in which most scans take it is not the search the user asked
        /// for. Reporting the count lets that be noticed rather than inferred from a disappointing result.
        /// </para>
        /// </summary>
        /// <param name="scans">The scans about to be searched.</param>
        /// <param name="commonParameters">Parameters for this file; the match mode is read from here.</param>
        /// <param name="dataFileName">Optional file name, included in the message when supplied.</param>
        public static string GetMonoisotopicFallbackWarning(this IEnumerable<Ms2ScanWithSpecificMass> scans,
            CommonParameters commonParameters, string dataFileName = null)
        {
            if (!UsesMostAbundantPeak(commonParameters) || scans == null)
            {
                return null;
            }

            int total = 0;
            int withoutApex = 0;
            foreach (var scan in scans.Where(s => s != null))
            {
                total++;
                if (scan.PrecursorMostAbundantMass <= 0)
                {
                    withoutApex++;
                }
            }

            if (total == 0 || withoutApex == 0)
            {
                return null;
            }

            string where = string.IsNullOrWhiteSpace(dataFileName) ? "" : $" in {dataFileName}";
            double percent = 100.0 * withoutApex / total;
            return $"Most-abundant precursor matching: {withoutApex} of {total} MS2 scans ({percent:F1}%){where} "
                 + "have no observed most-abundant peak and will be matched on their monoisotopic precursor "
                 + "mass instead. This occurs when no isotopic envelope was deconvoluted for the precursor — "
                 + "for example a scan-header precursor, or a neutral mass read from a pre-deconvoluted file.";
        }

        public static bool UsesMostAbundantPeak(this CommonParameters commonParameters)
            => commonParameters?.PrecursorMassMatchMode == PrecursorMassMatchMode.MostAbundant;

        public static double ApexOffset(this CommonParameters commonParameters, double monoisotopicMass)
        {
            AverageResidue averagine = commonParameters.GetAverageResidue();
            return averagine.GetDiffToMonoisotopic(averagine.GetMostIntenseMassIndex(monoisotopicMass));
        }

        public static double ApexOffset(this AverageResidue residue, double monoisotopicMass)
            => residue.GetDiffToMonoisotopic(residue.GetMostIntenseMassIndex(monoisotopicMass));

        public static AverageResidue GetAverageResidue(this CommonParameters commonParameters)
            => commonParameters.PrecursorDeconvolutionParameters?.AverageResidueModel ?? DefaultAveragine;

        public static double IsotopeSpacing(this CommonParameters commonParameters)
            => commonParameters.PrecursorDeconvolutionParameters?.ExpectedIsotopeSpacing ?? Constants.C13MinusC12;
    }
}
