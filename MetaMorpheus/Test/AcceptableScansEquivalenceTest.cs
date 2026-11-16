using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Util;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    /// <summary>
    /// The Classic search finds the scans in a peptide's precursor windows with ClassicSearchEngine.FindAcceptableScans: a
    /// bucket table (SortedMassIndex) in place of a binary search over every scan's mass, the intervals of a plain
    /// DotMassDiffAcceptor read without allocating them, and the scans collected into a list rather than yielded. These
    /// tests hold it to the previous iterator, kept verbatim below as of 38423d28c: the same (notch, scan) pairs in the same
    /// order, for every kind of search mode.
    /// </summary>
    [TestFixture]
    internal static class AcceptableScansEquivalenceTest
    {
        [Test]
        public static void FirstAtOrAboveIsWhatBinarySearchFound()
        {
            var random = new Random(77);
            long compared = 0;
            for (int trial = 0; trial < 3000; trial++)
            {
                double[] masses = RandomAscendingMasses(random, trial);
                var index = new SortedMassIndex(masses);
                foreach (double value in Queries(random, masses))
                {
                    Assert.That(index.FirstAtOrAbove(value), Is.EqualTo(GetFirstScanWithMassOverOrEqualBefore(masses, value)),
                        $"trial {trial}, value {value:R}");
                    compared++;
                }
            }
            Assert.That(compared, Is.GreaterThan(500_000));
        }

        [Test]
        public static void FindAcceptableScansGivesTheScansOfTheOldIterator()
        {
            var random = new Random(1618);
            MsDataScan ms2 = new TestDataFile().GetAllScansList().First(s => s.MsnOrder == 2);
            long hits = 0;

            for (int trial = 0; trial < 300; trial++)
            {
                var commonParameters = new CommonParameters();
                Ms2ScanWithSpecificMass[] scans = RandomScans(random, ms2, commonParameters, 1 + random.Next(trial % 10 == 0 ? 5 : 3000));
                double[] scanMasses = scans.Select(b => b.GetPrecursorMassForSearch(commonParameters)).ToArray(); // as the engine builds it
                MassDiffAcceptor searchMode = RandomSearchMode(random, trial);
                var engine = new ClassicSearchEngine(new SpectralMatch[scans.Length], scans, new List<Modification>(), new List<Modification>(),
                    null, null, null, new List<IBioPolymer>(), searchMode, commonParameters, null, null, new List<string>(), false);
                var actual = new List<ScanWithIndexAndNotchInfo>();

                for (int p = 0; p < 200; p++)
                {
                    // peptide masses anywhere, and right on top of scans so that windows hold several
                    double peptideMass = random.Next(3) == 0
                        ? random.NextDouble() * 6000
                        : scanMasses[random.Next(scanMasses.Length)] + (random.NextDouble() - 0.5) * 1e-4;

                    var expected = GetAcceptableScansBefore(peptideMass, searchMode, scans, scanMasses).Select(s => (s.Notch, s.ScanIndex)).ToList();
                    engine.FindAcceptableScans(peptideMass, actual);
                    Assert.That(actual.Select(s => (s.Notch, s.ScanIndex)), Is.EqualTo(expected), $"trial {trial} ({searchMode.GetType().Name}), peptide {peptideMass:R}");
                    hits += expected.Count;
                }
            }
            Assert.That(hits, Is.GreaterThan(50_000));
        }

        private static MassDiffAcceptor RandomSearchMode(Random random, int trial)
        {
            var tolerance = new PpmTolerance(1 + random.Next(20));
            // notches close enough together that windows overlap, so a scan can be listed under two notches
            double[] shifts = Enumerable.Range(0, 1 + random.Next(80)).Select(_ => Math.Round((random.NextDouble() - 0.3) * 200, random.Next(2) == 0 ? 2 : 6))
                .Append(0).Distinct().ToArray();
            return (trial % 6) switch
            {
                0 or 1 => new DotMassDiffAcceptor("", shifts, tolerance),
                2 => new SinglePpmAroundZeroSearchMode(1 + random.Next(20)),
                3 => new SingleAbsoluteAroundZeroSearchMode(random.NextDouble()),
                4 => new IntervalMassDiffAcceptor("", shifts.Select(s => new DoubleRange(s - 0.01, s + random.NextDouble() * 0.5))),
                _ => trial % 12 == 5 ? new OpenSearchMode() : new MostAbundantDotMassDiffAcceptor("", shifts, tolerance, new Averagine(), 2),
            };
        }

        /// <summary> Scans sorted by search mass as the tasks sort them, with repeated precursors (equal masses). </summary>
        private static Ms2ScanWithSpecificMass[] RandomScans(Random random, MsDataScan ms2, CommonParameters commonParameters, int count)
        {
            var scans = new List<Ms2ScanWithSpecificMass>();
            while (scans.Count < count)
            {
                int charge = 1 + random.Next(6);
                double mz = random.Next(20) == 0 ? 2000 + random.NextDouble() * 4000 : 300 + random.NextDouble() * 1200;
                int repeats = random.Next(10) == 0 ? 2 + random.Next(3) : 1;
                for (int r = 0; r < repeats; r++)
                {
                    scans.Add(new Ms2ScanWithSpecificMass(ms2, mz, charge, "file", commonParameters, new IsotopicEnvelope[0]));
                }
            }
            return scans.OrderBy(b => b.GetPrecursorMassForSearch(commonParameters)).ToArray();
        }

        /// <summary> Ascending, as a sorted precursor list is: repeats, zeros, an occasional NaN-first, negative or infinite array. </summary>
        private static double[] RandomAscendingMasses(Random random, int trial)
        {
            int count = trial % 40 == 0 ? random.Next(2) : 1 + random.Next(trial % 7 == 0 ? 20 : 4000);
            double maxMass = trial % 9 == 0 ? 200_000 : 6000;
            var masses = new List<double>();
            while (masses.Count < count)
            {
                double mass = random.Next(10) == 0 ? Math.Round(random.NextDouble() * maxMass * 2) / 2 : random.NextDouble() * maxMass; // some on bucket edges
                int repeats = random.Next(10) == 0 ? 2 + random.Next(3) : 1;
                masses.AddRange(Enumerable.Repeat(mass, repeats));
            }
            if (trial % 25 == 3) masses.Add(double.NaN);
            if (trial % 25 == 4) masses.Add(-random.NextDouble() * 10);
            if (trial % 25 == 5) masses.Add(double.PositiveInfinity);
            if (trial % 25 == 6) masses.Add(0.0);
            if (trial % 25 == 7) masses.Add(-0.0);
            return masses.OrderBy(m => m).ToArray();
        }

        private static IEnumerable<double> Queries(Random random, double[] masses)
        {
            yield return double.NaN;
            yield return double.PositiveInfinity;
            yield return double.NegativeInfinity;
            yield return 0.0;
            yield return -0.0;
            yield return -5;
            for (int q = 0; q < 200; q++)
            {
                if (masses.Length == 0 || random.Next(5) == 0)
                {
                    yield return random.NextDouble() * 7000;
                    continue;
                }
                double m = masses[random.Next(masses.Length)];
                switch (random.Next(7))
                {
                    case 0: yield return m; break;
                    case 1: yield return Math.BitIncrement(m); break;
                    case 2: yield return Math.BitDecrement(m); break;
                    case 3: yield return Math.Floor(m * 2) / 2; break;                 // on a half-dalton bucket edge
                    case 4: yield return Math.Floor(m / 64) * 64; break;               // on a wide bucket edge
                    case 5: yield return m + (random.NextDouble() - 0.5) * 0.01; break;
                    default: yield return masses[^1] + random.NextDouble() * 100; break;
                }
            }
        }

        // ===== verbatim as of 38423d28c (ClassicSearchEngine.GetAcceptableScans and GetFirstScanWithMassOverOrEqual), with the
        // engine's fields passed in

        private sealed class ScanWithIndexAndNotchInfoBefore(int notch, int scanIndex)
        {
            public readonly int Notch = notch;
            public readonly int ScanIndex = scanIndex;
        }

        private static IEnumerable<ScanWithIndexAndNotchInfoBefore> GetAcceptableScansBefore(double peptideMonoisotopicMass, MassDiffAcceptor searchMode,
            Ms2ScanWithSpecificMass[] ArrayOfSortedMS2Scans, double[] MyScanPrecursorMasses)
        {
            foreach (AllowedIntervalWithNotch allowedIntervalWithNotch in searchMode.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(peptideMonoisotopicMass))
            {
                int scanIndex = GetFirstScanWithMassOverOrEqualBefore(MyScanPrecursorMasses, allowedIntervalWithNotch.Minimum);
                if (scanIndex < ArrayOfSortedMS2Scans.Length)
                {
                    var scanMass = MyScanPrecursorMasses[scanIndex];
                    while (scanMass <= allowedIntervalWithNotch.Maximum)
                    {
                        yield return new ScanWithIndexAndNotchInfoBefore(allowedIntervalWithNotch.Notch, scanIndex);
                        scanIndex++;
                        if (scanIndex == ArrayOfSortedMS2Scans.Length)
                        {
                            break;
                        }

                        scanMass = MyScanPrecursorMasses[scanIndex];
                    }
                }
            }
        }

        private static int GetFirstScanWithMassOverOrEqualBefore(double[] MyScanPrecursorMasses, double minimum)
        {
            int index = Array.BinarySearch(MyScanPrecursorMasses, minimum);
            if (index < 0)
            {
                index = ~index;
            }

            // index of the first element that is larger than value
            return index;
        }
    }
}
