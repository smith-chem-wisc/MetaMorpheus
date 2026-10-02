using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace EngineLayer.Util
{
    /// <summary>
    /// Orders peptide ids by their coarse (indexed) score, highest first, optionally keeping only the best.
    ///
    /// Coarse scores are bytes, so a counting sort over 256 buckets replaces a comparison sort: two linear
    /// passes, no per-scan allocation. It is stable, deliberately: peptides that tie on the coarse score are
    /// handed back in the order they were observed, which is what OrderByDescending did. Ties decide which
    /// peptide ends up as the PSM when the fine scores also tie, so an unstable sort would quietly change results.
    ///
    /// This is the only place the search sorts or cuts candidates by coarse score; modern, crosslink and glyco
    /// all come through here.
    ///
    /// One instance per thread -- it reuses its buffers and is not safe to share.
    /// </summary>
    public sealed class DescendingScoreSorter
    {
        private readonly int[] _countsByScore = new int[byte.MaxValue + 1];
        private readonly List<int> _sorted = new List<int>();

        /// <summary>
        /// Returns <paramref name="peptideIds"/> ordered by descending score, ties in their existing order.
        /// </summary>
        /// <remarks>
        /// The span points into a buffer this instance reuses, so the next call on the same instance overwrites
        /// it. Finish iterating before sorting again; a nested sort needs its own sorter.
        /// </remarks>
        public ReadOnlySpan<int> Sort(List<int> peptideIds, ScanScoringTable scores)
        {
            SelectTop(peptideIds, scores, 0, 0, _sorted);
            return CollectionsMarshal.AsSpan(_sorted);
        }

        /// <summary>
        /// Keeps the candidate peptides worth matching: every id scoring at least <paramref name="scoreCutoff"/>, cut after the
        /// <paramref name="topN"/>th best but keeping all ids tied with it, written to <paramref name="topCandidates"/> from highest score
        /// to lowest and, within a score, in the order the ids were observed. That is exactly what a stable descending sort by score,
        /// stopped below the topN-th score, produces. A <paramref name="topN"/> of zero or less keeps every id at or above the cutoff.
        /// </summary>
        /// <param name="topCandidates"> Cleared and filled. </param>
        public void SelectTop(List<int> candidateIds, ScanScoringTable scores, int scoreCutoff, int topN, List<int> topCandidates)
        {
            topCandidates.Clear();
            foreach (int id in candidateIds)
            {
                _countsByScore[scores[id]]++;
            }

            // The lowest score kept: the score of the topN-th best candidate, or the cutoff when there are fewer than topN.
            int lowestKeptScore = Math.Max(scoreCutoff, 0);
            if (topN > 0)
            {
                int atOrAbove = 0;
                for (int score = byte.MaxValue; score >= lowestKeptScore; score--)
                {
                    atOrAbove += _countsByScore[score];
                    if (atOrAbove >= topN)
                    {
                        lowestKeptScore = score;
                        break;
                    }
                }
            }

            // Where each kept score's ids start in the output, highest score first.
            int kept = 0;
            for (int score = byte.MaxValue; score >= lowestKeptScore; score--)
            {
                int count = _countsByScore[score];
                _countsByScore[score] = kept;
                kept += count;
            }

            CollectionsMarshal.SetCount(topCandidates, kept);
            Span<int> output = CollectionsMarshal.AsSpan(topCandidates);
            foreach (int id in candidateIds)
            {
                int score = scores[id];
                if (score >= lowestKeptScore)
                {
                    output[_countsByScore[score]++] = id;
                }
            }

            Array.Clear(_countsByScore);
        }
    }
}
