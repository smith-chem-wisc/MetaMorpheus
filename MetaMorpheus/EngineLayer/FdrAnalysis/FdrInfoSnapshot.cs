using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer.SpectrumMatch;

namespace EngineLayer.FdrAnalysis
{
    /// <summary>
    /// Shields the FDR state of a set of spectral matches from an <see cref="FdrAnalysisEngine"/> run over a
    /// subset of them, and puts it back on <see cref="Dispose"/>.
    ///
    /// The engine writes its results onto the matches it is given: PsmFdrInfo and PeptideFdrInfo, the
    /// notch q-values and cumulative counts on each hypothesis, and PsmCount. A second pass over part of
    /// the matches (the individual-file results, which re-run FDR file by file) would otherwise leave its
    /// values behind for every later reader that expects the first pass's.
    ///
    /// Taking the snapshot gives each match a copy of its FdrInfo objects, so the subset pass writes into
    /// the copies, starting from exactly the values it would have seen before (including PEP, which it
    /// reads but does not recompute). The original FdrInfo objects are never written to.
    /// </summary>
    public sealed class FdrInfoSnapshot : IDisposable
    {
        private readonly List<(SpectralMatch Match, FdrInfo Psm, FdrInfo Peptide, int PsmCount)> _matches = new();
        private readonly List<(SpectralMatchHypothesis Hypothesis, double? QValueNotch, double? CumulativeTargetNotch, double? CumulativeDecoyNotch,
            double? PeptideQValueNotch, double? PeptideCumulativeTargetNotch, double? PeptideCumulativeDecoyNotch)> _hypotheses = new();
        private bool _restored;

        private FdrInfoSnapshot(IEnumerable<SpectralMatch> matches)
        {
            foreach (var match in matches.Where(m => m != null))
            {
                _matches.Add((match, match.PsmFdrInfo, match.PeptideFdrInfo, match.PsmCount));
                match.PsmFdrInfo = Copy(match.PsmFdrInfo);
                match.PeptideFdrInfo = Copy(match.PeptideFdrInfo);

                foreach (var h in match.BestMatchingBioPolymersWithSetMods)
                {
                    _hypotheses.Add((h, h.QValueNotch, h.CumulativeTargetNotch, h.CumulativeDecoyNotch,
                        h.PeptideQValueNotch, h.PeptideCumulativeTargetNotch, h.PeptideCumulativeDecoyNotch));
                }
            }
        }

        /// <summary>Snapshots the FDR state of <paramref name="matches"/>; dispose the result to restore it.</summary>
        public static FdrInfoSnapshot Take(IEnumerable<SpectralMatch> matches)
        {
            ArgumentNullException.ThrowIfNull(matches);
            return new FdrInfoSnapshot(matches);
        }

        /// <summary>Puts back the FDR state every match had when the snapshot was taken.</summary>
        public void Dispose()
        {
            if (_restored) return;
            _restored = true;

            foreach (var (match, psm, peptide, psmCount) in _matches)
            {
                match.PsmFdrInfo = psm;
                match.PeptideFdrInfo = peptide;
                match.PsmCount = psmCount;
            }

            foreach (var s in _hypotheses)
            {
                s.Hypothesis.QValueNotch = s.QValueNotch;
                s.Hypothesis.CumulativeTargetNotch = s.CumulativeTargetNotch;
                s.Hypothesis.CumulativeDecoyNotch = s.CumulativeDecoyNotch;
                s.Hypothesis.PeptideQValueNotch = s.PeptideQValueNotch;
                s.Hypothesis.PeptideCumulativeTargetNotch = s.PeptideCumulativeTargetNotch;
                s.Hypothesis.PeptideCumulativeDecoyNotch = s.PeptideCumulativeDecoyNotch;
            }
        }

        private static FdrInfo Copy(FdrInfo fdrInfo) => fdrInfo == null ? null : new FdrInfo
        {
            CumulativeTarget = fdrInfo.CumulativeTarget,
            CumulativeDecoy = fdrInfo.CumulativeDecoy,
            CumulativeTargetNotch = fdrInfo.CumulativeTargetNotch,
            CumulativeDecoyNotch = fdrInfo.CumulativeDecoyNotch,
            QValue = fdrInfo.QValue,
            QValueNotch = fdrInfo.QValueNotch,
            PEP = fdrInfo.PEP,
            PEP_QValue = fdrInfo.PEP_QValue
        };
    }
}
