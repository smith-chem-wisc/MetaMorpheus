using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    /// <summary>
    /// The Classic search scores a candidate with MetaMorpheusEngine.ScoreFragmentIons, which builds the matched-ion list
    /// only when the score reaches the cutoff. These tests hold it to what the search did before - MatchFragmentIons and
    /// then CalculatePeptideScore, both kept verbatim below as of 38423d28c - bit for bit: the same score always, and the
    /// same ions, in the same order, from the same objects, whenever the score reaches the cutoff.
    /// </summary>
    [TestFixture]
    internal static class ScoreFragmentIonsEquivalenceTest
    {
        [Test]
        public static void ScoreFragmentIonsGivesTheScoreAndIonsOfMatchingThenScoring()
        {
            var random = new Random(3141);
            List<MsDataScan> ms2Scans = new TestDataFile().GetAllScansList().Where(s => s.MsnOrder == 2).ToList();
            var scratch = new FragmentMatchScratch();
            int listsCompared = 0;
            int nullsCompared = 0;
            int nonZeroScores = 0;

            for (int trial = 0; trial < 4000; trial++)
            {
                Tolerance tolerance = trial % 2 == 0
                    ? new PpmTolerance(1 + random.Next(30))
                    : new AbsoluteTolerance(0.001 + random.NextDouble() * 0.05);
                // LowCID leaves the scan without deconvoluted fragments, the other path that matches nothing
                var dissociation = trial % 40 == 0 ? DissociationType.LowCID : DissociationType.HCD;
                var commonParameters = new CommonParameters(productMassTolerance: tolerance, dissociationType: dissociation);
                int precursorCharge = (1 + random.Next(5)) * (random.Next(8) == 0 ? -1 : 1);
                MsDataScan ms2 = ms2Scans[random.Next(ms2Scans.Count)];
                IsotopicEnvelope[] envelopes = RandomEnvelopes(random, trial % 25 == 0 ? 0 : 1 + random.Next(200));
                var scan = new Ms2ScanWithSpecificMass(ms2, 1000.0.ToMz(Math.Abs(precursorCharge)), precursorCharge, "file",
                    commonParameters, envelopes);
                List<Product> products = RandomProducts(random, envelopes, tolerance, 1 + random.Next(80));
                Assert.That(MetaMorpheusEngine.CanScoreBeforeMatching(scan, commonParameters, matchAllCharges: false), Is.True);

                List<MatchedFragmentIon> expectedIons = MatchFragmentIonsBefore(scan, products, commonParameters);
                double expectedScore = CalculatePeptideScoreBefore(scan.TheScan, expectedIons);
                if (expectedScore != 0)
                {
                    nonZeroScores++;
                }

                // cutoffs below, at and just above the score, and anywhere in the range a search uses
                double cutoff = (trial % 4) switch
                {
                    0 => expectedScore,
                    1 => Math.BitIncrement(expectedScore),
                    2 => random.NextDouble() * 4,
                    _ => -1 + random.NextDouble() * 12
                };

                double actualScore = MetaMorpheusEngine.ScoreFragmentIons(scan, products, commonParameters, cutoff, scratch, out var actualIons);

                Assert.That(BitConverter.DoubleToInt64Bits(actualScore), Is.EqualTo(BitConverter.DoubleToInt64Bits(expectedScore)), $"trial {trial}: score");
                if (expectedScore >= cutoff)
                {
                    Assert.That(actualIons, Is.Not.Null, $"trial {trial}: ions at or above the cutoff");
                    AssertIdentical(expectedIons, actualIons, trial);
                    listsCompared++;
                }
                else
                {
                    Assert.That(actualIons, Is.Null, $"trial {trial}: no ions below the cutoff");
                    nullsCompared++;
                }
            }

            Assert.That(listsCompared, Is.GreaterThan(1000));
            Assert.That(nullsCompared, Is.GreaterThan(1000));
            Assert.That(nonZeroScores, Is.GreaterThan(3000));
        }

        [Test]
        public static void ScoringBeforeMatchingIsOnlyOfferedOnThePlainPath()
        {
            MsDataScan ms2 = new TestDataFile().GetAllScansList().First(s => s.MsnOrder == 2);
            var plain = new CommonParameters();
            var scan = new Ms2ScanWithSpecificMass(ms2, 500, 2, "file", plain);

            Assert.That(MetaMorpheusEngine.CanScoreBeforeMatching(scan, plain, matchAllCharges: false), Is.True);
            Assert.That(MetaMorpheusEngine.CanScoreBeforeMatching(scan, plain, matchAllCharges: true), Is.False);
            Assert.That(MetaMorpheusEngine.CanScoreBeforeMatching(scan, new CommonParameters(addCompIons: true), matchAllCharges: false), Is.False);
        }

        [Test]
        public static void ScratchGrowsForLongProductLists()
        {
            var random = new Random(5);
            MsDataScan ms2 = new TestDataFile().GetAllScansList().First(s => s.MsnOrder == 2);
            var commonParameters = new CommonParameters();
            IsotopicEnvelope[] envelopes = RandomEnvelopes(random, 300);
            var scan = new Ms2ScanWithSpecificMass(ms2, 1000.0.ToMz(3), 3, "file", commonParameters, envelopes);
            // every product lands on an envelope, so all of them are matched and the scratch must hold them all
            var products = envelopes.Select((e, i) => new Product(ProductType.b, FragmentationTerminus.N, e.MonoisotopicMass, i + 1, i + 1, 0)).ToList();

            var expectedIons = MatchFragmentIonsBefore(scan, products, commonParameters);
            double score = MetaMorpheusEngine.ScoreFragmentIons(scan, products, commonParameters, 0, new FragmentMatchScratch(), out var ions);

            Assert.That(BitConverter.DoubleToInt64Bits(score), Is.EqualTo(BitConverter.DoubleToInt64Bits(CalculatePeptideScoreBefore(scan.TheScan, expectedIons))));
            AssertIdentical(expectedIons, ions, 0);
            Assert.That(ions.Count, Is.GreaterThan(64));
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
                    var peaks = new List<(double mz, double intensity)> { (mass.ToMz(Math.Abs(charge)), intensity) };
                    if (random.Next(3) == 0)
                    {
                        peaks.Add(((mass + 1.00335).ToMz(Math.Abs(charge)), intensity * random.NextDouble()));
                    }
                    envelopes.Add(new IsotopicEnvelope(peaks, mass, charge, intensity, 0));
                }
            }
            return envelopes.OrderBy(p => p.MonoisotopicMass).ToArray();
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

        // ===== verbatim as of 38423d28c. The closest-envelope search is the old one, run on the mass array exactly as the
        // old constructor built it, so nothing below goes through code this branch has changed.

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

        /// <summary> MetaMorpheusEngine.MatchFragmentIons, matchAllCharges false, no envelopes, AddCompIons false. </summary>
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

        /// <summary> MetaMorpheusEngine.CalculatePeptideScore, fragmentsCanHaveDifferentCharges false. </summary>
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
