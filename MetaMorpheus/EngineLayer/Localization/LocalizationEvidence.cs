using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Localization
{
    /// <summary>
    /// A modification placed on one residue of a peptide.
    ///
    /// A mass on a peptide terminus is in exactly the fragments that hold the terminal residue, so it fragments as the
    /// same mass on that residue does, and no spectrum can tell the two apart. The first residue's site therefore stands
    /// for that residue or the peptide N-terminus, and the last residue's site for that residue or the C-terminus.
    /// </summary>
    /// <param name="ResidueIndex">Zero-based index of the residue in the peptide.</param>
    /// <param name="Residue">The residue at <paramref name="ResidueIndex"/>.</param>
    /// <param name="Modification">The modification placed there, or null for an unidentified mass difference.</param>
    public readonly record struct LocalizationSite(int ResidueIndex, char Residue, Modification Modification = null);

    /// <summary>
    /// One candidate arrangement of a match's modifications: the site of each modification whose position is in
    /// question. Modifications whose position is not in question (fixed modifications, terminal modifications) are not
    /// listed.
    /// </summary>
    public sealed class LocalizationPlacement
    {
        public LocalizationPlacement(IReadOnlyList<LocalizationSite> sites)
        {
            Sites = sites;
        }

        /// <summary>The sites, in residue order.</summary>
        public IReadOnlyList<LocalizationSite> Sites { get; }

        public override string ToString() => string.Join(",", Sites.Select(s => $"{s.Residue}{s.ResidueIndex + 1}"));
    }

    /// <summary>
    /// The spectral evidence for each candidate placement of a PSM's modifications, collected while the spectra file
    /// is loaded so that localization probabilities can be computed later, after every file has been read.
    ///
    /// Holds every distinct experimental peak matched by at least one candidate, as its intensity relative to the
    /// spectrum's base peak, and for each candidate which of those peaks it matched. A peak matched by every candidate
    /// says nothing about the site, but is kept: the peaks of unmodified PSMs are what a model of matched-peak
    /// intensities is fitted from.
    /// </summary>
    public sealed class LocalizationEvidence
    {
        private readonly float[] _peakRelativeIntensities;
        private readonly ulong[] _matchMasks; // candidate c's bits are words [c * _wordsPerCandidate, (c + 1) * _wordsPerCandidate)
        private readonly int _wordsPerCandidate;
        private readonly int[] _predictedFragmentCounts;

        internal LocalizationEvidence(double massDifference, int peptideLength, bool isProteinNTerminal, bool isProteinCTerminal,
            LocalizationPlacement[] placements, double[] scores, float[] peakRelativeIntensities, ulong[] matchMasks,
            int[] predictedFragmentCounts = null, float[] noiseCoverage = null, int sharedPredictedFragmentCount = -1)
        {
            SharedPredictedFragmentCount = sharedPredictedFragmentCount;
            MassDifference = massDifference;
            PeptideLength = peptideLength;
            IsProteinNTerminal = isProteinNTerminal;
            IsProteinCTerminal = isProteinCTerminal;
            Placements = placements;
            Scores = scores;
            _peakRelativeIntensities = peakRelativeIntensities;
            _matchMasks = matchMasks;
            _wordsPerCandidate = WordsFor(peakRelativeIntensities.Length);
            _predictedFragmentCounts = predictedFragmentCounts;
            NoiseCoverage = noiseCoverage;
        }

        /// <summary>
        /// The observed minus the theoretical mass of the match's peptide. For placements of an unidentified mass
        /// difference this is the mass placed; for placements of identified modifications it is near zero.
        /// </summary>
        public double MassDifference { get; }

        /// <summary>The number of residues in the peptide.</summary>
        public int PeptideLength { get; }

        /// <summary>Whether the peptide begins its protein (or begins it after the initiator methionine is removed).</summary>
        public bool IsProteinNTerminal { get; }

        /// <summary>Whether the peptide ends its protein.</summary>
        public bool IsProteinCTerminal { get; }

        /// <summary>The candidate placements, in the order their evidence is indexed.</summary>
        public IReadOnlyList<LocalizationPlacement> Placements { get; }

        /// <summary>The MetaMorpheus score of the peptide with each candidate placement.</summary>
        public IReadOnlyList<double> Scores { get; }

        public int CandidateCount => Placements.Count;

        /// <summary>The number of distinct experimental peaks matched by at least one candidate.</summary>
        public int PeakCount => _peakRelativeIntensities.Length;

        /// <summary>The intensity of a matched peak divided by the intensity of the spectrum's base peak.</summary>
        public double PeakRelativeIntensity(int peak) => _peakRelativeIntensities[peak];

        /// <summary>Whether the candidate at <paramref name="candidate"/> matched the peak at <paramref name="peak"/>.</summary>
        public bool Matches(int candidate, int peak)
        {
            ulong word = _matchMasks[candidate * _wordsPerCandidate + (peak >> 6)];
            return (word & (1UL << (peak & 63))) != 0;
        }

        /// <summary>Whether every candidate matched the peak, so that it cannot tell the candidates apart.</summary>
        public bool IsSharedByAllCandidates(int peak)
        {
            for (int c = 0; c < CandidateCount; c++)
            {
                if (!Matches(c, peak))
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>Whether the number of fragments each candidate predicted within the scan's range was recorded.</summary>
        public bool HasPredictedFragmentCounts => _predictedFragmentCounts != null;

        /// <summary>
        /// The number of theoretical fragments the candidate predicts within the scan's range, matched or not, for the
        /// missing-fragment term of <see cref="LocalizationLikelihood.IntensityRatio"/>.
        /// </summary>
        public int PredictedFragmentCount(int candidate) => _predictedFragmentCounts[candidate];

        /// <summary>
        /// The number of fragments within the scan's range that every candidate predicts alike (the same neutral mass),
        /// so that they are predicted whichever candidate is right; -1 when not recorded. With the peaks every candidate
        /// matches, it lets a modified match teach how often a predicted fragment is observed without knowing its site.
        /// </summary>
        public int SharedPredictedFragmentCount { get; }

        /// <summary>
        /// The scan's analyzer and activation, e.g. "IonTrap2D CID", or null when not recorded. One file can hold scans of
        /// several types (a Fusion method alternating ion-trap and Orbitrap MS2), whose fragment statistics differ.
        /// </summary>
        public string ScanType { get; init; }

        /// <summary>
        /// The peaks no candidate explains, as a histogram over <see cref="PeakIntensityRatioModel"/>'s intensity bins
        /// in which each peak adds the fraction of the fragment mass range its tolerance window covers: bin b holds the
        /// chance that a predicted fragment lands on an unexplained peak of that intensity. Null when not recorded.
        /// </summary>
        public IReadOnlyList<float> NoiseCoverage { get; }

        internal static int WordsFor(int peakCount) => (peakCount + 63) >> 6;
    }

    /// <summary>
    /// Collects the matched fragment ions of each candidate placement of one PSM and builds its
    /// <see cref="LocalizationEvidence"/>. One builder serves one PSM on one thread.
    /// </summary>
    internal sealed class LocalizationEvidenceBuilder
    {
        private readonly double _basePeakIntensity;
        private readonly Dictionary<(double Mz, double Intensity), int> _peakIndex = new();
        private readonly List<float> _peakRelativeIntensities = new();
        private readonly List<LocalizationPlacement> _placements = new();
        private readonly List<double> _scores = new();
        private readonly List<List<int>> _matchedPeaks = new();
        private readonly List<int> _predictedFragmentCounts = new();
        private float[] _noiseCoverage;
        private int _sharedPredictedFragments = -1;
        private readonly string _scanType;

        /// <param name="basePeakIntensity">The intensity of the most intense peak in the spectrum.</param>
        /// <param name="scanType">The scan's analyzer and activation (see <see cref="LocalizationEvidence.ScanType"/>).</param>
        /// <param name="perCleavage">Whether a candidate's ions are kept one per backbone cleavage: of a fragment and its
        /// neutral-loss forms (and its charge states), only the most intense matched peak.</param>
        public LocalizationEvidenceBuilder(double basePeakIntensity, string scanType = null, bool perCleavage = false)
        {
            _basePeakIntensity = basePeakIntensity;
            _scanType = scanType;
            _perCleavage = perCleavage;
        }

        private readonly bool _perCleavage;
        private readonly HashSet<(double Mz, double Intensity)> _explained = new();

        /// <summary>Of each cleavage's matched ions, the most intense: one piece of evidence per cleavage.</summary>
        private static IEnumerable<MatchedFragmentIon> OnePerCleavage(IReadOnlyList<MatchedFragmentIon> ions) => ions
            .GroupBy(ion => (ion.NeutralTheoreticalProduct.ProductType,
                (long)Math.Round((ion.NeutralTheoreticalProduct.NeutralMass + ion.NeutralTheoreticalProduct.NeutralLoss) * 1e4)))
            .Select(g => g.MaxBy(ion => ion.Intensity));

        /// <summary>
        /// Records a candidate placement and the ions it matched. Ions with the same m/z and intensity are one
        /// experimental peak (two products can match one isotopic envelope), so a peak is counted once per candidate.
        /// </summary>
        /// <param name="predictedFragments">The number of theoretical fragments the candidate predicts within the scan's
        /// range, or -1 when not known; the counts are kept only when every candidate gives one.</param>
        public void AddCandidate(LocalizationPlacement placement, double score, IReadOnlyList<MatchedFragmentIon> matchedIons, int predictedFragments = -1)
        {
            var peaks = new List<int>(matchedIons.Count);
            foreach (MatchedFragmentIon ion in matchedIons)
            {
                _explained.Add((Math.Round(ion.Mz, 6), ion.Intensity)); // explained, even where not kept as evidence
            }
            foreach (MatchedFragmentIon ion in _perCleavage ? OnePerCleavage(matchedIons) : matchedIons)
            {
                var key = (Math.Round(ion.Mz, 6), ion.Intensity);
                if (!_peakIndex.TryGetValue(key, out int index))
                {
                    index = _peakRelativeIntensities.Count;
                    _peakIndex.Add(key, index);
                    _peakRelativeIntensities.Add((float)(ion.Intensity / _basePeakIntensity));
                }
                peaks.Add(index);
            }

            _placements.Add(placement);
            _scores.Add(score);
            _matchedPeaks.Add(peaks);
            _predictedFragmentCounts.Add(predictedFragments);
        }

        /// <summary>Sets <see cref="LocalizationEvidence.SharedPredictedFragmentCount"/>.</summary>
        public void SetSharedPredictedFragments(int count) => _sharedPredictedFragments = count;

        /// <summary>Sets the unexplained-peak coverage histogram directly (see <see cref="LocalizationEvidence.NoiseCoverage"/>).</summary>
        public void SetNoiseCoverage(float[] coverage)
        {
            if (coverage.Length != PeakIntensityRatioModel.BinCount)
            {
                throw new ArgumentException($"Expected {PeakIntensityRatioModel.BinCount} bins.", nameof(coverage));
            }
            _noiseCoverage = coverage;
        }

        /// <summary>
        /// Records, once every candidate is added, the isotopic envelopes of the scan that no candidate matched: each
        /// adds to its intensity bin the width of its tolerance window over the span of fragment masses that can be
        /// observed, so that bin b holds the chance a predicted fragment lands on an unexplained envelope of that
        /// intensity. Envelopes outside the span, or of a charge the matching does not admit, are left out.
        /// </summary>
        public void RecordUnmatchedEnvelopes(IReadOnlyList<IsotopicEnvelope> envelopes, Tolerance tolerance, double lowestMass, double highestMass, int maxCharge)
        {
            var coverage = new float[PeakIntensityRatioModel.BinCount];
            double span = highestMass - lowestMass;
            if (envelopes != null && span > 0)
            {
                foreach (IsotopicEnvelope envelope in envelopes)
                {
                    double mass = envelope.MonoisotopicMass;
                    if (mass < lowestMass || mass > highestMass || Math.Abs(envelope.Charge) > Math.Abs(maxCharge))
                    {
                        continue;
                    }
                    // the key a matched envelope has (see MetaMorpheusEngine.MatchFragmentIons)
                    double intensity = envelope.Peaks.First().intensity;
                    if (_explained.Contains((Math.Round(mass.ToMz(envelope.Charge), 6), intensity)))
                    {
                        continue;
                    }
                    double width = tolerance.GetMaximumValue(mass) - tolerance.GetMinimumValue(mass);
                    coverage[PeakIntensityRatioModel.BinOf(intensity / _basePeakIntensity)] += (float)(width / span);
                }
            }
            _noiseCoverage = coverage;
        }

        public LocalizationEvidence Build(double massDifference, int peptideLength, bool isProteinNTerminal, bool isProteinCTerminal)
        {
            int words = LocalizationEvidence.WordsFor(_peakRelativeIntensities.Count);
            var masks = new ulong[words * _placements.Count];
            for (int c = 0; c < _matchedPeaks.Count; c++)
            {
                foreach (int peak in _matchedPeaks[c])
                {
                    masks[c * words + (peak >> 6)] |= 1UL << (peak & 63);
                }
            }

            int[] predicted = _predictedFragmentCounts.All(n => n >= 0) ? _predictedFragmentCounts.ToArray() : null;
            return new LocalizationEvidence(massDifference, peptideLength, isProteinNTerminal, isProteinCTerminal,
                _placements.ToArray(), _scores.ToArray(), _peakRelativeIntensities.ToArray(), masks, predicted, _noiseCoverage, _sharedPredictedFragments) { ScanType = _scanType };
        }
    }
}
