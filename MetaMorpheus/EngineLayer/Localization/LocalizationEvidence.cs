using Omics.Fragmentation;
using System;
using System.Collections.Generic;

namespace EngineLayer.Localization
{
    /// <summary>
    /// One candidate placement of a modification on a peptide: a residue.
    ///
    /// A mass on a peptide terminus is in exactly the fragments that hold the terminal residue, so it fragments as the
    /// same mass on that residue does, and no spectrum can tell the two apart. The first residue's site therefore stands
    /// for that residue or the peptide N-terminus, and the last residue's site for that residue or the C-terminus.
    /// </summary>
    /// <param name="ResidueIndex">Zero-based index of the residue in the peptide.</param>
    /// <param name="Residue">The amino acid at <paramref name="ResidueIndex"/>.</param>
    public readonly record struct LocalizationSite(int ResidueIndex, char Residue);

    /// <summary>
    /// The spectral evidence for each candidate placement of a PSM's modification, collected while the spectra file is
    /// loaded so that localization probabilities can be computed later, after every file has been read.
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

        internal LocalizationEvidence(double massDifference, bool isProteinNTerminal, bool isProteinCTerminal,
            LocalizationSite[] sites, double[] scores, float[] peakRelativeIntensities, ulong[] matchMasks)
        {
            MassDifference = massDifference;
            IsProteinNTerminal = isProteinNTerminal;
            IsProteinCTerminal = isProteinCTerminal;
            Sites = sites;
            Scores = scores;
            _peakRelativeIntensities = peakRelativeIntensities;
            _matchMasks = matchMasks;
            _wordsPerCandidate = WordsFor(peakRelativeIntensities.Length);
        }

        /// <summary>The mass placed at each candidate site (observed minus theoretical peptide mass).</summary>
        public double MassDifference { get; }

        /// <summary>Whether the peptide begins its protein (or begins it after the initiator methionine is removed).</summary>
        public bool IsProteinNTerminal { get; }

        /// <summary>Whether the peptide ends its protein.</summary>
        public bool IsProteinCTerminal { get; }

        /// <summary>The candidate placements, in the order their evidence is indexed.</summary>
        public IReadOnlyList<LocalizationSite> Sites { get; }

        /// <summary>The MetaMorpheus score of the peptide with the modification at each candidate site.</summary>
        public IReadOnlyList<double> Scores { get; }

        public int CandidateCount => Sites.Count;

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
        private readonly List<LocalizationSite> _sites = new();
        private readonly List<double> _scores = new();
        private readonly List<List<int>> _matchedPeaks = new();

        /// <param name="basePeakIntensity">The intensity of the most intense peak in the spectrum.</param>
        public LocalizationEvidenceBuilder(double basePeakIntensity)
        {
            _basePeakIntensity = basePeakIntensity;
        }

        /// <summary>
        /// Records a candidate placement and the ions it matched. Ions with the same m/z and intensity are one
        /// experimental peak (two products can match one isotopic envelope), so a peak is counted once per candidate.
        /// </summary>
        public void AddCandidate(LocalizationSite site, double score, IReadOnlyList<MatchedFragmentIon> matchedIons)
        {
            var peaks = new List<int>(matchedIons.Count);
            foreach (MatchedFragmentIon ion in matchedIons)
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
            _sites.Add(site);
            _scores.Add(score);
            _matchedPeaks.Add(peaks);
        }

        public LocalizationEvidence Build(double massDifference, bool isProteinNTerminal, bool isProteinCTerminal)
        {
            int words = LocalizationEvidence.WordsFor(_peakRelativeIntensities.Count);
            var masks = new ulong[words * _sites.Count];
            for (int c = 0; c < _matchedPeaks.Count; c++)
            {
                foreach (int peak in _matchedPeaks[c])
                {
                    masks[c * words + (peak >> 6)] |= 1UL << (peak & 63);
                }
            }

            return new LocalizationEvidence(massDifference, isProteinNTerminal, isProteinCTerminal,
                _sites.ToArray(), _scores.ToArray(), _peakRelativeIntensities.ToArray(), masks);
        }
    }
}
