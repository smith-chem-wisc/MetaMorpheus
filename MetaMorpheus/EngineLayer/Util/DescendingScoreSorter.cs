using System;
using System.Collections.Generic;

namespace EngineLayer.Util
{
    /// <summary>
    /// Orders peptide ids by their coarse (indexed) score, highest first.
    ///
    /// The search used to do this with OrderByDescending, once per scan, which allocated a LINQ
    /// ordering plus a key array every time. Coarse scores are bytes, so a counting sort over 256
    /// buckets does the same job in two linear passes with no per-scan allocation.
    ///
    /// It is stable, deliberately: peptides that tie on the coarse score are handed back in the order
    /// they were found, which is what OrderByDescending did. Ties decide which peptide ends up as the
    /// PSM when the fine scores also tie, so an unstable sort here would quietly change results.
    ///
    /// One instance per thread -- it reuses its buffers and is not safe to share.
    /// </summary>
    public sealed class DescendingScoreSorter
    {
        private const int BucketCount = byte.MaxValue + 1;

        private readonly int[] _bucketStarts = new int[BucketCount + 1];
        private int[] _sorted = new int[64];

        /// <summary>
        /// Returns <paramref name="peptideIds"/> ordered by descending score in <paramref name="scores"/>,
        /// ties in their existing order. The span is valid until the next call on this instance.
        /// </summary>
        public ReadOnlySpan<int> Sort(List<int> peptideIds, ScanScoringTable scores)
        {
            int count = peptideIds.Count;
            if (count == 0)
            {
                return ReadOnlySpan<int>.Empty;
            }

            if (_sorted.Length < count)
            {
                _sorted = new int[Math.Max(count, _sorted.Length * 2)];
            }

            Array.Clear(_bucketStarts, 0, _bucketStarts.Length);

            // Bucket by inverted score, so bucket 0 holds the highest scores.
            for (int i = 0; i < count; i++)
            {
                _bucketStarts[byte.MaxValue - scores[peptideIds[i]] + 1]++;
            }

            for (int bucket = 0; bucket < BucketCount; bucket++)
            {
                _bucketStarts[bucket + 1] += _bucketStarts[bucket];
            }

            // Walking the ids in their original order keeps equal scores in that order.
            for (int i = 0; i < count; i++)
            {
                int id = peptideIds[i];
                _sorted[_bucketStarts[byte.MaxValue - scores[id]]++] = id;
            }

            return new ReadOnlySpan<int>(_sorted, 0, count);
        }
    }
}
