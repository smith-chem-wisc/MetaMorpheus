using EngineLayer.GlycoSearch;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    /// <summary>
    /// Glyco search sends only each scan's TopN coarse candidates on to glycan matching. When the database is
    /// split into partitions that cut has to be taken over the whole database, not per partition, or the
    /// identifications depend on the partition count -- which RaisePartitionsToFitMemory can change on its own.
    /// </summary>
    [TestFixture]
    public static class GlycoPartitionInvarianceTest
    {
        /// <summary>
        /// Pooling each partition's own TopN cut and cutting again must select exactly the candidates a single
        /// partition would. Scores are drawn from a narrow range so ties at the cut are common, since keeping every
        /// tie with the TopN-th score is where a pooled cut could most easily diverge.
        /// </summary>
        [Test]
        public static void PooledPartitionCuts_SelectTheSameCandidatesAsOnePartition()
        {
            var random = new Random(20260915);
            int[] topNs = { 0, 1, 3, 5, 50 };
            int trialsWithMoreCandidatesThanTopN = 0;

            for (int trial = 0; trial < 2000; trial++)
            {
                int peptideCount = random.Next(1, 300);
                byte cutoff = (byte)random.Next(0, 4);
                int topN = topNs[random.Next(topNs.Length)];
                int partitions = random.Next(2, 9);

                var scores = new byte[peptideCount];
                for (int i = 0; i < peptideCount; i++)
                {
                    scores[i] = (byte)random.Next(0, 12);
                }

                // observed order is arbitrary in the engine, so shuffle it rather than use id order
                List<int> observed = Enumerable.Range(0, peptideCount).OrderBy(_ => random.Next()).ToList();

                var singlePartition = new List<int>();
                GlycoSearchEngine.SelectTopN(observed, scores, cutoff, topN, singlePartition);
                if (topN > 0 && observed.Count(id => scores[id] >= cutoff) > topN)
                {
                    trialsWithMoreCandidatesThanTopN++;
                }

                // contiguous partitions, as the task slices the protein list
                var pooled = new List<(int Partition, int PeptideId, byte Score)>();
                for (int partition = 0; partition < partitions; partition++)
                {
                    int start = partition * peptideCount / partitions;
                    int end = (partition + 1) * peptideCount / partitions;
                    List<int> observedHere = observed.Where(id => id >= start && id < end).ToList();

                    var cutHere = new List<int>();
                    GlycoSearchEngine.SelectTopN(observedHere, scores, cutoff, topN, cutHere);
                    if (cutHere.Count == 0)
                    {
                        continue;
                    }

                    pooled.AddRange(cutHere.Select(id => (partition, id, scores[id])));
                    pooled = GlycoSearchEngine.KeepGlobalTopN(pooled, topN);
                }

                Assert.That(pooled.Select(c => c.PeptideId).OrderBy(id => id), Is.EqualTo(singlePartition.OrderBy(id => id)),
                    $"trial {trial}: {peptideCount} peptides, cutoff {cutoff}, TopN {topN}, {partitions} partitions");
                Assert.That(pooled.Select(c => c.Score), Is.Ordered.Descending, $"trial {trial}: candidates are not highest score first");
            }

            Assert.That(trialsWithMoreCandidatesThanTopN, Is.GreaterThan(500), "premise: the cut must actually remove candidates in most trials");
        }

        [Test]
        public static void KeepGlobalTopN_KeepsEveryTieWithTheTopNthScore()
        {
            var candidates = new List<(int Partition, int PeptideId, byte Score)>
            {
                (0, 0, 5), (0, 1, 9), (1, 0, 7), (1, 1, 7), (2, 0, 7), (2, 1, 3),
            };

            var kept = GlycoSearchEngine.KeepGlobalTopN(candidates, 2);

            Assert.That(kept, Is.EqualTo(new List<(int, int, byte)> { (0, 1, 9), (1, 0, 7), (1, 1, 7), (2, 0, 7) }),
                "the second-highest score is 7, so all three 7s stay, in the order they were added");
        }

        [Test]
        public static void KeepGlobalTopN_NonPositiveTopNKeepsEverything()
        {
            var candidates = new List<(int Partition, int PeptideId, byte Score)> { (0, 0, 1), (1, 0, 4), (1, 1, 2) };

            Assert.That(GlycoSearchEngine.KeepGlobalTopN(candidates, 0).Count, Is.EqualTo(3));
            Assert.That(GlycoSearchEngine.KeepGlobalTopN(candidates, -1).Count, Is.EqualTo(3));
        }
    }
}
