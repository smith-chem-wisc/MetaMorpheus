using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    /// <summary>
    /// The Classic search sorts a peptide's product masses once (SortedProductMasses) and finds the closest envelope of
    /// every product in one forward walk of each scan's mass array (Ms2ScanWithSpecificMass.GetClosestFragmentMasses)
    /// instead of one binary search per product. These tests hold the walk to the previous binary search, and the scoring
    /// built on it to the previous MatchFragmentIons + CalculatePeptideScore, all kept verbatim below as of 38423d28c.
    /// </summary>
    [TestFixture]
    internal static class MergedWalkEquivalenceTest
    {
        [Test]
        public static void WalkFindsTheClosestEnvelopeTheBinarySearchFound()
        {
            var random = new Random(99);
            MsDataScan ms2 = new TestDataFile().GetAllScansList().First(s => s.MsnOrder == 2);
            var commonParameters = new CommonParameters();
            long compared = 0;

            for (int trial = 0; trial < 3000; trial++)
            {
                IsotopicEnvelope[] envelopes = RandomEnvelopesWithOddities(random, trial);
                var scan = new Ms2ScanWithSpecificMass(ms2, 500, 2, "file", commonParameters, envelopes);
                double[] masses = envelopes.Select(p => p.MonoisotopicMass).ToArray(); // as the constructor builds it

                double[] queries = Queries(random, masses, trial);
                var closest = new int[queries.Length];
                scan.GetClosestFragmentMasses(queries, queries.Length, closest);

                for (int k = 0; k < queries.Length; k++)
                {
                    int expected = masses.Length == 0 ? -1 : GetClosestFragmentMassBefore(masses, queries[k]);
                    Assert.That(closest[k], Is.EqualTo(expected), $"trial {trial}, query {k} = {queries[k]:R}");
                    compared++;
                }
            }
            Assert.That(compared, Is.GreaterThan(300_000));
        }

        [Test]
        public static void SortedProductMassesAreTheProductsInMassOrder()
        {
            var random = new Random(11);
            var sorted = new SortedProductMasses();
            for (int trial = 0; trial < 500; trial++)
            {
                var products = Enumerable.Range(0, random.Next(150))
                    .Select(i => new Product(ProductType.b, FragmentationTerminus.N,
                        random.Next(10) == 0 ? double.NaN : Math.Round(random.NextDouble() * 2000, 1), i + 1, i + 1, 0))
                    .ToList();
                sorted.Build(products);

                var expected = products.Where(p => !double.IsNaN(p.NeutralMass)).Select(p => p.NeutralMass).OrderBy(m => m).ToArray();
                Assert.That(sorted.Count, Is.EqualTo(expected.Length));
                Assert.That(sorted.ProductCount, Is.EqualTo(products.Count));
                Assert.That(sorted.Masses.Take(sorted.Count), Is.EqualTo(expected));
                Assert.That(sorted.ProductIndices.Take(sorted.Count).OrderBy(i => i),
                    Is.EqualTo(Enumerable.Range(0, products.Count).Where(i => !double.IsNaN(products[i].NeutralMass))));
                for (int k = 0; k < sorted.Count; k++)
                {
                    Assert.That(products[sorted.ProductIndices[k]].NeutralMass, Is.EqualTo(sorted.Masses[k]));
                }
            }
        }

        [Test]
        public static void ScoringWithTheWalkGivesTheScoreAndIonsOfMatchingThenScoring()
        {
            var random = new Random(2718);
            List<MsDataScan> ms2Scans = new TestDataFile().GetAllScansList().Where(s => s.MsnOrder == 2).ToList();
            var scratch = new FragmentMatchScratch();
            var sorted = new SortedProductMasses();
            int listsCompared = 0;

            for (int trial = 0; trial < 4000; trial++)
            {
                Tolerance tolerance = trial % 2 == 0
                    ? new PpmTolerance(1 + random.Next(30))
                    : new AbsoluteTolerance(0.001 + random.NextDouble() * 0.05);
                var dissociation = trial % 40 == 0 ? DissociationType.LowCID : DissociationType.HCD;
                var commonParameters = new CommonParameters(productMassTolerance: tolerance, dissociationType: dissociation);
                int precursorCharge = (1 + random.Next(5)) * (random.Next(8) == 0 ? -1 : 1);
                IsotopicEnvelope[] envelopes = RandomEnvelopes(random, trial % 25 == 0 ? 0 : 1 + random.Next(200));
                var scan = new Ms2ScanWithSpecificMass(ms2Scans[random.Next(ms2Scans.Count)], 1000.0.ToMz(Math.Abs(precursorCharge)),
                    precursorCharge, "file", commonParameters, envelopes);

                // real peptides' products (b and y interleaved, as Fragment returns them) or random ones with NaN and diagnostic ions
                List<Product> products = trial % 3 == 0 ? RealProducts(random) : RandomProducts(random, envelopes, tolerance, 1 + random.Next(80));
                if (trial % 3 == 0 && envelopes.Length > 0)
                {
                    // put some envelopes on the real products so that they match
                    envelopes = envelopes.Concat(products.Where(p => random.Next(3) == 0).Select(p => Envelope(random, p.NeutralMass)))
                        .OrderBy(e => e.MonoisotopicMass).ToArray();
                    scan = new Ms2ScanWithSpecificMass(scan.TheScan, scan.PrecursorMonoisotopicPeakMz, precursorCharge, "file", commonParameters, envelopes);
                }
                sorted.Build(products);

                List<MatchedFragmentIon> expectedIons = MatchFragmentIonsBefore(scan, products, commonParameters);
                double expectedScore = CalculatePeptideScoreBefore(scan.TheScan, expectedIons);
                double cutoff = (trial % 4) switch
                {
                    0 => expectedScore,
                    1 => Math.BitIncrement(expectedScore),
                    _ => -1 + random.NextDouble() * 12
                };

                double actualScore = MetaMorpheusEngine.ScoreFragmentIons(scan, products, sorted, commonParameters, cutoff, scratch, out var actualIons);

                Assert.That(BitConverter.DoubleToInt64Bits(actualScore), Is.EqualTo(BitConverter.DoubleToInt64Bits(expectedScore)), $"trial {trial}: score");
                if (expectedScore >= cutoff)
                {
                    AssertIdentical(expectedIons, actualIons, trial);
                    listsCompared++;
                }
                else
                {
                    Assert.That(actualIons, Is.Null, $"trial {trial}: no ions below the cutoff");
                }
            }
            Assert.That(listsCompared, Is.GreaterThan(1500));
        }

        [Test]
        public static void SortedMassesOfAnotherProductListAreNotUsed()
        {
            var random = new Random(4);
            MsDataScan ms2 = new TestDataFile().GetAllScansList().First(s => s.MsnOrder == 2);
            var commonParameters = new CommonParameters();
            IsotopicEnvelope[] envelopes = RandomEnvelopes(random, 150);
            var scan = new Ms2ScanWithSpecificMass(ms2, 1000.0.ToMz(2), 2, "file", commonParameters, envelopes);
            List<Product> products = RandomProducts(random, envelopes, commonParameters.ProductMassTolerance, 40);
            var stale = new SortedProductMasses();
            stale.Build(products.Take(10).ToList());

            var expectedIons = MatchFragmentIonsBefore(scan, products, commonParameters);
            double score = MetaMorpheusEngine.ScoreFragmentIons(scan, products, stale, commonParameters, 0, new FragmentMatchScratch(), out var ions);

            Assert.That(BitConverter.DoubleToInt64Bits(score), Is.EqualTo(BitConverter.DoubleToInt64Bits(CalculatePeptideScoreBefore(scan.TheScan, expectedIons))));
            AssertIdentical(expectedIons, ions, 0);
        }

        private static List<Product> RealProducts(Random random)
        {
            string[] sequences = { "PEPTIDEK", "LGVEGQKHR", "SAMPLERMAK", "ACDEFGHIKLMNPQRSTVWY", "GGGGGGK", "MKAAQDNSQR" };
            var peptide = new PeptideWithSetModifications(sequences[random.Next(sequences.Length)], new Dictionary<string, Omics.Modifications.Modification>());
            var products = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            return products;
        }

        private static IsotopicEnvelope Envelope(Random random, double mass)
        {
            int charge = 1 + random.Next(3);
            double intensity = 1 + random.NextDouble() * 1e6;
            return new IsotopicEnvelope(new List<(double mz, double intensity)> { (mass.ToMz(charge), intensity) }, mass, charge, intensity, 0);
        }

        private static IsotopicEnvelope[] RandomEnvelopes(Random random, int count)
        {
            var envelopes = new List<IsotopicEnvelope>();
            for (int i = 0; i < count; i++)
            {
                double mass = 100 + random.NextDouble() * 3500;
                int repeats = random.Next(10) == 0 ? 2 + random.Next(3) : 1;
                for (int r = 0; r < repeats; r++)
                {
                    int charge = (1 + random.Next(8)) * (random.Next(10) == 0 ? -1 : 1);
                    double intensity = 1 + random.NextDouble() * 1e6;
                    envelopes.Add(new IsotopicEnvelope(new List<(double mz, double intensity)> { (mass.ToMz(Math.Abs(charge)), intensity) },
                        mass, charge, intensity, 0));
                }
            }
            return envelopes.OrderBy(p => p.MonoisotopicMass).ToArray();
        }

        /// <summary> As in ClosestFragmentMassEquivalenceTest: repeats, NaN, signed zeros, infinities, some out of order. </summary>
        private static IsotopicEnvelope[] RandomEnvelopesWithOddities(Random random, int trial)
        {
            int count = trial % 50 == 0 ? 0 : 1 + random.Next(250);
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
                    < 20 => Math.Round(100 + random.NextDouble() * 3500, 1),
                    _ => 100 + random.NextDouble() * 3500
                };
                if (trial % 3 != 0 && (double.IsNaN(mass) || double.IsInfinity(mass) || mass <= 0))
                {
                    mass = 100 + random.NextDouble() * 3500;
                }
                int repeats = random.Next(8) == 0 ? 2 + random.Next(3) : 1;
                for (int r = 0; r < repeats; r++)
                {
                    envelopes.Add(new IsotopicEnvelope(new List<(double mz, double intensity)> { (500, 1000) }, mass, 1, 1000, 0));
                }
            }
            return trial % 11 == 5
                ? envelopes.OrderBy(_ => random.Next()).ToArray()
                : envelopes.OrderBy(p => p.MonoisotopicMass).ToArray();
        }

        /// <summary>
        /// Ascending query masses as a sorted product list gives them - with repeats, exact hits, ties between neighbours and
        /// the ends - and, now and then, NaN or masses out of order, which the walk must still answer correctly.
        /// </summary>
        private static double[] Queries(Random random, double[] masses, int trial)
        {
            var queries = new List<double>();
            for (int q = 0; q < 120; q++)
            {
                if (masses.Length == 0 || random.Next(6) == 0)
                {
                    queries.Add(random.NextDouble() * 4000 - 100);
                    continue;
                }
                int k = random.Next(masses.Length);
                double m = masses[k];
                switch (random.Next(7))
                {
                    case 0: queries.Add(m); break;
                    case 1: queries.Add(Math.BitIncrement(m)); break;
                    case 2: queries.Add(Math.BitDecrement(m)); break;
                    case 3: queries.Add(k > 0 ? (m + masses[k - 1]) / 2 : m - 1); break;
                    case 4: queries.Add(m + (random.NextDouble() - 0.5) * 0.02); break;
                    case 5: queries.Add(m); queries.Add(m); break;
                    default: queries.Add(random.NextDouble() * 4000); break;
                }
            }
            queries.Add(double.PositiveInfinity);
            queries.Add(double.NegativeInfinity);
            var ordered = queries.OrderBy(x => x).ToList();
            if (trial % 13 == 0)
            {
                ordered.Insert(random.Next(ordered.Count), double.NaN);
            }
            if (trial % 17 == 0)
            {
                ordered = ordered.OrderBy(_ => random.Next()).ToList();
            }
            return ordered.ToArray();
        }

        private static List<Product> RandomProducts(Random random, IsotopicEnvelope[] envelopes, Tolerance tolerance, int count)
        {
            ProductType[] types = { ProductType.b, ProductType.y, ProductType.c, ProductType.zDot, ProductType.D, ProductType.M };
            var products = new List<Product>();
            for (int i = 0; i < count; i++)
            {
                double mass;
                int kind = random.Next(10);
                if (envelopes.Length == 0 || kind < 2)
                {
                    mass = 100 + random.NextDouble() * 3500;
                }
                else
                {
                    double target = envelopes[random.Next(envelopes.Length)].MonoisotopicMass;
                    double halfWidth = (tolerance.GetMaximumValue(target) - tolerance.GetMinimumValue(target)) / 2;
                    mass = kind switch
                    {
                        2 or 3 => target,
                        4 or 5 => target + halfWidth * (random.NextDouble() * 1.8 - 0.9),
                        6 or 7 => target + halfWidth * (1.01 + random.NextDouble()) * (random.Next(2) == 0 ? -1 : 1),
                        8 => target + halfWidth * (random.Next(2) == 0 ? -1 : 1),
                        _ => double.NaN
                    };
                }
                products.Add(new Product(types[random.Next(types.Length)], FragmentationTerminus.N, mass, 1 + i, 1 + i, 0));
            }
            return products;
        }

        private static void AssertIdentical(List<MatchedFragmentIon> expected, List<MatchedFragmentIon> actual, int trial)
        {
            Assert.That(actual, Is.Not.Null, $"trial {trial}: ions");
            Assert.That(actual.Count, Is.EqualTo(expected.Count), $"trial {trial}: ion count");
            for (int i = 0; i < expected.Count; i++)
            {
                var e = expected[i];
                var a = actual[i];
                Assert.That(a.GetType(), Is.EqualTo(e.GetType()), $"trial {trial}, ion {i}: type");
                Assert.That(a.NeutralTheoreticalProduct, Is.SameAs(e.NeutralTheoreticalProduct), $"trial {trial}, ion {i}: product");
                Assert.That(BitConverter.DoubleToInt64Bits(a.Mz), Is.EqualTo(BitConverter.DoubleToInt64Bits(e.Mz)), $"trial {trial}, ion {i}: m/z");
                Assert.That(BitConverter.DoubleToInt64Bits(a.Intensity), Is.EqualTo(BitConverter.DoubleToInt64Bits(e.Intensity)), $"trial {trial}, ion {i}: intensity");
                Assert.That(a.Charge, Is.EqualTo(e.Charge), $"trial {trial}, ion {i}: charge");
            }
        }

        // ===== verbatim as of 38423d28c, closest-envelope search included, run on the mass array as the old constructor built it

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

        private static IsotopicEnvelope GetClosestExperimentalIsotopicEnvelopeBefore(Ms2ScanWithSpecificMass scan, double theoreticalNeutralMass)
        {
            double[] DeconvolutedMonoisotopicMasses = scan.ExperimentalFragments != null && scan.ExperimentalFragments.Any()
                ? scan.ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray()
                : new double[0];
            if (DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return null;
            }
            return scan.ExperimentalFragments[GetClosestFragmentMassBefore(DeconvolutedMonoisotopicMasses, theoreticalNeutralMass)];
        }

        private static List<MatchedFragmentIon> MatchFragmentIonsBefore(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, CommonParameters commonParameters)
        {
            var productMassTolerance = commonParameters.ProductMassTolerance;
            var matchedFragmentIons = new List<MatchedFragmentIon>();

            if (scan.TheScan.MassSpectrum.XcorrProcessed && scan.TheScan.MassSpectrum.XArray.Length != 0)
            {
                throw new NotSupportedException("not the path under test");
            }

            // if the spectrum has no peaks
            if (scan.ExperimentalFragments != null && !scan.ExperimentalFragments.Any())
            {
                return matchedFragmentIons;
            }

            // search for ions in the spectrum
            for (int i = 0; i < theoreticalProducts.Count; i++)
            {
                var product = theoreticalProducts[i];
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(product.NeutralMass))
                {
                    continue;
                }

                // get the closest peak in the spectrum to the theoretical peak
                var closestExperimentalMass = GetClosestExperimentalIsotopicEnvelopeBefore(scan, product.NeutralMass);

                // is the mass error acceptable?
                if (closestExperimentalMass != null
                    && productMassTolerance.Within(closestExperimentalMass.MonoisotopicMass, product.NeutralMass)
                    && Math.Abs(closestExperimentalMass.Charge) <= Math.Abs(scan.PrecursorCharge))
                {
                    matchedFragmentIons.Add(new MatchedFragmentIon(product, closestExperimentalMass.MonoisotopicMass.ToMz(closestExperimentalMass.Charge),
                        closestExperimentalMass.Peaks.First().intensity, closestExperimentalMass.Charge));
                }
            }

            return matchedFragmentIons;
        }

        private static double CalculatePeptideScoreBefore(MsDataScan thisScan, List<MatchedFragmentIon> matchedFragmentIons)
        {
            double score = 0;

            if (thisScan.MassSpectrum.XcorrProcessed)
            {
                // XCorr
                foreach (var fragment in matchedFragmentIons)
                {
                    switch (fragment.NeutralTheoreticalProduct.ProductType)
                    {
                        case ProductType.aDegree:
                        case ProductType.aStar:
                        case ProductType.bWaterLoss:
                        case ProductType.bAmmoniaLoss:
                        case ProductType.yWaterLoss:
                        case ProductType.yAmmoniaLoss:
                            score += 0.01 * fragment.Intensity;
                            break;
                        case ProductType.D: //count nothing for diagnostic ions.
                            break;
                        default:
                            score += 1 * fragment.Intensity;
                            break;
                    }
                }
            }
            else
            {
                // Morpheus score
                for (int i = 0; i < matchedFragmentIons.Count; i++)
                {
                    switch (matchedFragmentIons[i].NeutralTheoreticalProduct.ProductType)
                    {
                        case ProductType.D:
                            break;
                        default:
                            score += 1 + matchedFragmentIons[i].Intensity / thisScan.TotalIonCurrent;
                            break;
                    }

                }
            }

            return score;
        }
    }
}
