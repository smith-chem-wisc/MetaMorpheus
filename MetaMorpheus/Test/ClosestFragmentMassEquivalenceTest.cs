using Chemistry;
using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    /// <summary>
    /// Ms2ScanWithSpecificMass.GetClosestFragmentMass finds the insertion point with an inlined lower-bound search instead
    /// of Array.BinarySearch, and hands exact hits, NaN masses and arrays that are not ascending numbers back to
    /// Array.BinarySearch. These tests hold it to the previous implementation, kept verbatim below, for every kind of
    /// query and array that can reach it.
    /// </summary>
    [TestFixture]
    internal static class ClosestFragmentMassEquivalenceTest
    {
        [Test]
        public static void ClosestFragmentMassIsUnchangedBySearchingWithoutBinarySearch()
        {
            var random = new Random(2024);
            MsDataScan ms2 = new TestDataFile().GetAllScansList().First(s => s.MsnOrder == 2);
            var commonParameters = new CommonParameters();
            long queries = 0;
            int fastArrays = 0;

            for (int trial = 0; trial < 3000; trial++)
            {
                IsotopicEnvelope[] envelopes = RandomEnvelopes(random, trial);
                var scan = new Ms2ScanWithSpecificMass(ms2, 500, 2, "file", commonParameters, envelopes);

                // exactly how the constructor builds the array the search runs on
                double[] masses = envelopes.Select(p => p.MonoisotopicMass).ToArray();
                if (masses.Length > 0 && !double.IsNaN(masses[0]) && masses.Zip(masses.Skip(1), (a, b) => b >= a).All(x => x))
                {
                    fastArrays++;
                }

                foreach (double query in Queries(random, masses))
                {
                    int expected = GetClosestFragmentMassBefore(masses, query);
                    int actual = scan.GetClosestFragmentMass(query);
                    Assert.That(actual, Is.EqualTo(expected), $"trial {trial}, query {query:R}");
                    queries++;

                    if (expected >= 0)
                    {
                        Assert.That(scan.GetClosestExperimentalIsotopicEnvelope(query), Is.SameAs(scan.ExperimentalFragments[expected]));
                    }
                }
            }

            Assert.That(queries, Is.GreaterThan(500_000));
            Assert.That(fastArrays, Is.GreaterThan(2000)); // most arrays take the new search, the rest exercise the fallback
        }

        /// <summary>
        /// Arrays as the constructor receives them: sorted by mass the way GetNeutralExperimentalFragments sorts them (so
        /// NaN first), with repeated masses, signed zeros and infinities; some empty; some deliberately out of order.
        /// </summary>
        private static IsotopicEnvelope[] RandomEnvelopes(Random random, int trial)
        {
            int count = trial % 50 == 0 ? 0 : trial % 7 == 0 ? 1 + random.Next(3) : 1 + random.Next(250);
            var envelopes = new List<IsotopicEnvelope>();
            for (int i = 0; i < count; i++)
            {
                double mass = random.Next(40) switch
                {
                    0 => double.NaN,
                    1 => double.PositiveInfinity,
                    2 => double.NegativeInfinity,
                    3 => 0.0,
                    4 => -0.0,
                    5 => -random.NextDouble() * 100,
                    // a coarse grid makes repeated masses common, as adjacent duplicates are in real spectra
                    < 20 => Math.Round(100 + random.NextDouble() * 3500, 1),
                    _ => 100 + random.NextDouble() * 3500
                };
                if (trial % 3 != 0 && (double.IsNaN(mass) || double.IsInfinity(mass) || mass <= 0))
                {
                    mass = 100 + random.NextDouble() * 3500; // most arrays hold only ordinary masses
                }
                int repeats = random.Next(8) == 0 ? 2 + random.Next(3) : 1;
                for (int r = 0; r < repeats; r++)
                {
                    int charge = 1 + random.Next(6);
                    envelopes.Add(new IsotopicEnvelope(new List<(double mz, double intensity)> { (500, 1000) }, mass, charge, 1000, 0));
                }
            }
            return trial % 11 == 5
                ? envelopes.OrderBy(_ => random.Next()).ToArray()
                : envelopes.OrderBy(p => p.MonoisotopicMass).ToArray();
        }

        private static IEnumerable<double> Queries(Random random, double[] masses)
        {
            yield return double.NaN;
            yield return double.PositiveInfinity;
            yield return double.NegativeInfinity;
            yield return 0.0;
            yield return -0.0;
            yield return -1;
            yield return 5000;
            for (int q = 0; q < 200; q++)
            {
                if (masses.Length == 0)
                {
                    yield return random.NextDouble() * 4000;
                    continue;
                }
                int k = random.Next(masses.Length);
                double m = masses[k];
                switch (random.Next(8))
                {
                    case 0: yield return m; break;                                                   // exact hit
                    case 1: yield return Math.BitIncrement(m); break;
                    case 2: yield return Math.BitDecrement(m); break;
                    case 3: yield return k > 0 ? (m + masses[k - 1]) / 2 : m - 1; break;             // a tie between neighbours
                    case 4: yield return m + (random.NextDouble() - 0.5) * 0.02; break;
                    case 5: yield return masses[0] - random.NextDouble() * 10; break;
                    case 6: yield return masses[^1] + random.NextDouble() * 10; break;
                    default: yield return random.NextDouble() * 4000; break;
                }
            }
        }

        /// <summary>
        /// Ms2ScanWithSpecificMass.GetClosestFragmentMass as of 38423d28c, verbatim, with the scan's mass array passed in.
        /// </summary>
        private static int GetClosestFragmentMassBefore(double[] DeconvolutedMonoisotopicMasses, double mass)
        {
            int index = Array.BinarySearch(DeconvolutedMonoisotopicMasses, mass);
            if (index >= 0)
            {
                return index;
            }
            index = ~index;

            if (index == DeconvolutedMonoisotopicMasses.Length)
            {
                return index - 1;
            }
            if (index == 0 || mass - DeconvolutedMonoisotopicMasses[index - 1] > DeconvolutedMonoisotopicMasses[index] - mass)
            {
                return index;
            }

            return index - 1;
        }
    }
}
