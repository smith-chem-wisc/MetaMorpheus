using Chemistry;
using EngineLayer;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.Util;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
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

        /// <summary>
        /// The window share is two binary searches per scan; this pins it to a linear count of the peptides in (min, max], the same
        /// collapsed window IndexScoreScan scores, for a multi-notch acceptor and a plain interval.
        /// </summary>
        [Test]
        public static void MeanWindowShareOfIndexMatchesALinearCount()
        {
            var peptideIndex = SortedTrypticIndex();
            var masses = peptideIndex.Where(p => !double.IsNaN(p.MonoisotopicMass)).Select(p => p.MonoisotopicMass).Where((_, i) => i % 7 == 0).ToList();
            var acceptors = new MassDiffAcceptor[]
            {
                SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), MassDiffAcceptorType.ThreeMM, null),
                new IntervalMassDiffAcceptor("wide", new[] { new DoubleRange(-50, 120) }),
            };

            foreach (MassDiffAcceptor acceptor in acceptors)
            {
                double expected = masses.Average(mass =>
                {
                    var notches = acceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(mass).ToList();
                    double low = notches.Min(n => n.Minimum);
                    double high = notches.Max(n => n.Maximum);
                    return peptideIndex.Count(p => p.MonoisotopicMass > low && p.MonoisotopicMass <= high) / (double)peptideIndex.Count;
                });

                double actual = ModernSearchEngine.MeanWindowShareOfIndex(peptideIndex, masses, acceptor);
                Assert.That(actual, Is.EqualTo(expected).Within(1e-12), acceptor.FileNameAddition);
            }

            Assert.That(ModernSearchEngine.MeanWindowShareOfIndex(peptideIndex, new double[0], acceptors[0]), Is.EqualTo(0));
        }

        /// <summary>
        /// Finite is not narrow. A wide Custom interval is bounded on both sides, so the acceptor check alone would stamp it, but its
        /// window holds far more of the index than stamping pays for; an ordinary tolerance holds far less.
        /// </summary>
        [Test]
        public static void AWideFiniteIntervalIsNotNarrowEnoughToStamp()
        {
            var peptideIndex = SortedTrypticIndex();
            var masses = peptideIndex.Where(p => !double.IsNaN(p.MonoisotopicMass)).Select(p => p.MonoisotopicMass).ToList();

            MassDiffAcceptor oneMissedMonoisotopic = SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), MassDiffAcceptorType.OneMM, null);
            MassDiffAcceptor wide = new IntervalMassDiffAcceptor("wide", new[] { new DoubleRange(-200, 500) });
            Assert.That(ScanScoringTable.IsWorthStamping(wide), Is.True, "bounded on both sides");

            Assert.That(ScanScoringTable.IsWindowNarrowEnoughToStamp(ModernSearchEngine.MeanWindowShareOfIndex(peptideIndex, masses, oneMissedMonoisotopic)), Is.True);
            Assert.That(ScanScoringTable.IsWindowNarrowEnoughToStamp(ModernSearchEngine.MeanWindowShareOfIndex(peptideIndex, masses, wide)), Is.False);
        }

        /// <summary>
        /// The engine's choice end to end, including the fallback: when the window share cannot be measured (no scans, or an index
        /// not sorted by mass) the acceptor alone decides, as it did before the share check. Open never stamps.
        /// </summary>
        [Test]
        public static void EngineStampsOnlyWhenTheWindowShareIsSmallOrCannotBeMeasured()
        {
            var sorted = SortedTrypticIndex();
            var unsorted = Enumerable.Reverse(sorted).ToList();
            var parameters = new CommonParameters();
            var scans = sorted.Where(p => !double.IsNaN(p.MonoisotopicMass)).Where((_, i) => i % 7 == 0)
                .Select(p => ScanAtMass(p.MonoisotopicMass, parameters)).ToArray();

            MassDiffAcceptor oneMissedMonoisotopic = SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), MassDiffAcceptorType.OneMM, null);
            MassDiffAcceptor wide = new IntervalMassDiffAcceptor("wide", new[] { new DoubleRange(-200, 500) });
            MassDiffAcceptor open = SearchTask.GetMassDiffAcceptor(new PpmTolerance(5), MassDiffAcceptorType.Open, null);

            Assert.That(new StampChoiceProbe(scans, sorted, oneMissedMonoisotopic, parameters).Stamps, Is.True, "narrow window");
            Assert.That(new StampChoiceProbe(scans, sorted, wide, parameters).Stamps, Is.False, "wide window, measured");
            Assert.That(new StampChoiceProbe(null, sorted, wide, parameters).Stamps, Is.True, "no scans: the acceptor decides");
            Assert.That(new StampChoiceProbe(scans, unsorted, wide, parameters).Stamps, Is.True, "unsorted index: the acceptor decides");
            Assert.That(new StampChoiceProbe(scans, sorted, open, parameters).Stamps, Is.False, "open never stamps");
        }

        private sealed class StampChoiceProbe : ModernSearchEngine
        {
            public StampChoiceProbe(Ms2ScanWithSpecificMass[] scans, List<IBioPolymerWithSetMods> peptideIndex, MassDiffAcceptor acceptor,
                CommonParameters parameters)
                : base(new SpectralMatch[scans?.Length ?? 0], scans, peptideIndex, null, 0, parameters,
                    new List<(string, CommonParameters)> { ("", parameters) }, acceptor, 0, new List<string>())
            {
            }

            public bool Stamps => UseStampedScoringTable;
        }

        private static Ms2ScanWithSpecificMass ScanAtMass(double mass, CommonParameters parameters)
        {
            var dataScan = new MsDataScan(new MzSpectrum(new[] { 500.0 }, new[] { 1000.0 }, false), 1, 2, true, Polarity.Positive, 1,
                new MzRange(0, 5000), "", MZAnalyzerType.Orbitrap, 1000, null, null, "");
            return new Ms2ScanWithSpecificMass(dataScan, mass.ToMz(1), 1, "", parameters, new IsotopicEnvelope[0]);
        }

        private static List<IBioPolymerWithSetMods> SortedTrypticIndex()
        {
            var random = new Random(20261006);
            const string residues = "ACDEFGHIKLMNPQRSTVWY";
            var digestion = new DigestionParams(protease: "trypsin", maxMissedCleavages: 2, minPeptideLength: 5);
            var peptides = Enumerable.Range(0, 400)
                .Select(i => new Protein(new string(Enumerable.Range(0, 300).Select(_ => residues[random.Next(residues.Length)]).ToArray()), "P" + i))
                .SelectMany(p => p.Digest(digestion, new List<Modification>(), new List<Modification>()))
                .Cast<IBioPolymerWithSetMods>()
                .ToList();
            return IndexingEngine.SortByMonoisotopicMass(peptides);
        }
    }
}
