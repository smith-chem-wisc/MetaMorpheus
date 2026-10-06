using EngineLayer;
using EngineLayer.Util;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// The per-scan coarse-score table and the sort over it. Nothing in the search tests pins tie order or exercises the stamp
    /// wrap, so a regression in either would still pass every search test while changing which peptide becomes the PSM.
    /// </summary>
    [TestFixture]
    public static class ScanScoringTableTests
    {
        /// <summary>
        /// Sort must hand back exactly what OrderByDescending did, ties included. Scores come from a small alphabet so that
        /// nearly every trial has ties, which is where an unstable sort would differ.
        /// </summary>
        [Test]
        public static void SortMatchesOrderByDescendingIncludingTies([Values(false, true)] bool stamped)
        {
            var random = new Random(20260831);
            var sorter = new DescendingScoreSorter();
            int trialsWithTies = 0;

            for (int trial = 0; trial < 2000; trial++)
            {
                int peptideCount = random.Next(1, 500);
                var table = new ScanScoringTable(peptideCount, stamped);
                table.BeginScan();

                int maxScore = trial % 4 == 0 ? 256 : random.Next(1, 8);
                for (int id = 0; id < peptideCount; id++)
                {
                    table.Set(id, (byte)random.Next(0, maxScore));
                }

                List<int> ids = Enumerable.Range(0, peptideCount).OrderBy(_ => random.Next()).Take(random.Next(0, peptideCount + 1)).ToList();
                if (ids.Select(id => table[id]).Distinct().Count() < ids.Count)
                {
                    trialsWithTies++;
                }

                int[] expected = ids.OrderByDescending(id => table[id]).ToArray();
                Assert.That(sorter.Sort(ids, table).ToArray(), Is.EqualTo(expected), $"trial {trial}");
            }

            Assert.That(trialsWithTies, Is.GreaterThan(1000), "premise: most trials must contain ties");
        }

        /// <summary>
        /// A stamped table must read exactly like a byte[] that is cleared every scan, across several wraps of the 255 stamps
        /// and with scores that themselves wrap past 255.
        /// </summary>
        [Test]
        public static void StampedTableReadsLikeAClearedByteArrayAcrossStampWraps()
        {
            var random = new Random(20260901);
            const int peptideCount = 64;
            var stamped = new ScanScoringTable(peptideCount, stamped: true);
            var reference = new byte[peptideCount];

            for (int scan = 0; scan < 255 * 4 + 7; scan++)
            {
                stamped.BeginScan();
                Array.Clear(reference);

                // a few peptides per scan, so most cells hold a stale stamp from an earlier scan
                for (int touch = 0; touch < 8; touch++)
                {
                    int id = random.Next(peptideCount);
                    if (random.Next(10) == 0)
                    {
                        byte score = (byte)random.Next(250, 256);
                        stamped.Set(id, score);
                        reference[id] = score;
                    }
                    else
                    {
                        Assert.That(stamped.Increment(id), Is.EqualTo(++reference[id]), $"scan {scan}, peptide {id}");
                    }
                }

                for (int id = 0; id < peptideCount; id++)
                {
                    Assert.That(stamped[id], Is.EqualTo(reference[id]), $"scan {scan}, peptide {id}");
                }
            }
        }

        /// <summary>
        /// Stamping is chosen only when the precursor window is bounded on both sides. ModOpen is the case that matters: it is
        /// not OpenSearchMode, but its [-187, +inf) interval leaves the window open below, so a scan touches every bin from its
        /// first peptide, and the GUI steers it towards modern search.
        /// </summary>
        [Test]
        [TestCase(MassDiffAcceptorType.Exact, true)]
        [TestCase(MassDiffAcceptorType.OneMM, true)]
        [TestCase(MassDiffAcceptorType.ThreeMM, true)]
        [TestCase(MassDiffAcceptorType.PlusOrMinusThreeMM, true)]
        [TestCase(MassDiffAcceptorType.MostAbundant_Exact, true)]
        [TestCase(MassDiffAcceptorType.MostAbundant_PlusMinusTwo, true)]
        [TestCase(MassDiffAcceptorType.ModOpen, false)]
        [TestCase(MassDiffAcceptorType.Open, false)]
        public static void StampingIsChosenOnlyForABoundedPrecursorWindow(MassDiffAcceptorType type, bool expected)
        {
            MassDiffAcceptor acceptor = SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), type, null);
            Assert.That(ScanScoringTable.IsWorthStamping(acceptor), Is.EqualTo(expected));
        }
    }
}
