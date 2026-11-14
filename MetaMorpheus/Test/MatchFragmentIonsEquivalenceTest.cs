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
    /// MatchFragmentIons tests a product's mass error on the scan's cached mass array before it reads the
    /// IsotopicEnvelope, because almost every product misses and dereferencing the envelope just to reject it
    /// costs a cache miss per product. These tests hold it to exactly what the previous implementation returned,
    /// which is kept verbatim below: the same ions, in the same order, built from the same objects, bit for bit.
    /// </summary>
    [TestFixture]
    internal static class MatchFragmentIonsEquivalenceTest
    {
        [Test]
        [TestCase(false, false, DissociationType.HCD)]
        [TestCase(true, false, DissociationType.HCD)]
        [TestCase(false, true, DissociationType.HCD)]
        [TestCase(true, true, DissociationType.HCD)]
        [TestCase(true, false, DissociationType.EThcD)]
        public static void MatchFragmentIonsIsUnchangedByCheckingTheMassFirst(bool addCompIons, bool includeEnvelope, DissociationType dissociationType)
        {
            var random = new Random(42);
            MsDataScan ms2 = new TestDataFile().GetAllScansList().First(s => s.MsnOrder == 2);
            int matched = 0;

            for (int trial = 0; trial < 400; trial++)
            {
                Tolerance tolerance = trial % 2 == 0
                    ? new PpmTolerance(1 + random.Next(30))
                    : new AbsoluteTolerance(0.001 + random.NextDouble() * 0.05);
                var commonParameters = new CommonParameters(productMassTolerance: tolerance, addCompIons: addCompIons,
                    dissociationType: dissociationType);
                int precursorCharge = (1 + random.Next(5)) * (random.Next(8) == 0 ? -1 : 1);
                IsotopicEnvelope[] envelopes = RandomEnvelopes(random, trial % 25 == 0 ? 0 : 1 + random.Next(200));
                var scan = new Ms2ScanWithSpecificMass(ms2, 1000.0.ToMz(Math.Abs(precursorCharge)), precursorCharge, "file",
                    commonParameters, envelopes);
                List<Product> products = RandomProducts(random, envelopes, tolerance, 1 + random.Next(60));

                List<MatchedFragmentIon> expected = MatchFragmentIonsBefore(scan, products, commonParameters, includeEnvelope);
                List<MatchedFragmentIon> actual = MetaMorpheusEngine.MatchFragmentIons(scan, products, commonParameters,
                    includeExperimentalEnvelope: includeEnvelope);

                AssertIdentical(expected, actual, trial);
                matched += actual.Count;
            }

            // the random products are built to land on, near and beside real envelopes, so most trials match something
            Assert.That(matched, Is.GreaterThan(1000));
        }

        [Test]
        public static void ClosestEnvelopeIndexAgreesWithClosestEnvelope()
        {
            var random = new Random(7);
            MsDataScan ms2 = new TestDataFile().GetAllScansList().First(s => s.MsnOrder == 2);
            var commonParameters = new CommonParameters();

            for (int trial = 0; trial < 200; trial++)
            {
                IsotopicEnvelope[] envelopes = RandomEnvelopes(random, trial % 20 == 0 ? 0 : 1 + random.Next(100));
                var scan = new Ms2ScanWithSpecificMass(ms2, 500, 2, "file", commonParameters, envelopes);

                for (int i = 0; i < 50; i++)
                {
                    double mass = random.NextDouble() * 4000;
                    if (envelopes.Length > 0 && i % 5 == 0)
                    {
                        mass = envelopes[random.Next(envelopes.Length)].MonoisotopicMass; // exact hits, including duplicated masses
                    }

                    IsotopicEnvelope closest = scan.GetClosestExperimentalIsotopicEnvelope(mass);
                    int index = scan.GetClosestExperimentalIsotopicEnvelopeIndex(mass);
                    if (closest == null)
                    {
                        Assert.That(index, Is.EqualTo(-1));
                        continue;
                    }
                    Assert.That(scan.ExperimentalFragments[index], Is.SameAs(closest));
                    Assert.That(BitConverter.DoubleToInt64Bits(scan.GetDeconvolutedMonoisotopicMass(index)),
                        Is.EqualTo(BitConverter.DoubleToInt64Bits(closest.MonoisotopicMass)));
                }
            }
        }

        /// <summary>
        /// Envelopes sorted by mass, as GetNeutralExperimentalFragments returns them, with some masses repeated
        /// (real spectra have adjacent duplicates, so which duplicate the binary search lands on matters) and a
        /// spread of charges on both sides of the precursor's.
        /// </summary>
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

        /// <summary>
        /// Products that hit an envelope exactly, sit just inside or just outside the tolerance of one, fall
        /// between envelopes, or have no mass at all; with diagnostic and precursor ions mixed in, which the
        /// complementary-ion loop skips.
        /// </summary>
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
                        4 or 5 => target + halfWidth * (random.NextDouble() * 1.8 - 0.9),        // inside
                        6 or 7 => target + halfWidth * (1.01 + random.NextDouble()) * (random.Next(2) == 0 ? -1 : 1), // outside
                        8 => target + halfWidth * (random.Next(2) == 0 ? -1 : 1),               // on the edge
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
                if (e is MatchedFragmentIonWithEnvelope withEnvelope)
                {
                    Assert.That(((MatchedFragmentIonWithEnvelope)a).Envelope, Is.SameAs(withEnvelope.Envelope), $"trial {trial}, ion {i}: envelope");
                }
            }
        }

        private static readonly Dictionary<DissociationType, List<double>> ComplementaryIonConversion = new()
        {
            { DissociationType.HCD, new List<double> { Constants.ProtonMass } },
            { DissociationType.EThcD, new List<double> { Constants.ProtonMass, 2 * Constants.ProtonMass } },
        };

        /// <summary>
        /// MetaMorpheusEngine.MatchFragmentIons as of 38423d28c, deconvoluted-fragment path, verbatim apart from the
        /// complementary-ion table, which is private to the engine and copied above for the two types tested.
        /// </summary>
        private static List<MatchedFragmentIon> MatchFragmentIonsBefore(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts,
            CommonParameters commonParameters, bool includeExperimentalEnvelope)
        {
            var productMassTolerance = commonParameters.ProductMassTolerance;
            var matchedFragmentIons = new List<MatchedFragmentIon>();

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
                var closestExperimentalMass = scan.GetClosestExperimentalIsotopicEnvelope(product.NeutralMass);

                // is the mass error acceptable?
                if (closestExperimentalMass != null
                    && productMassTolerance.Within(closestExperimentalMass.MonoisotopicMass, product.NeutralMass)
                    && Math.Abs(closestExperimentalMass.Charge) <= Math.Abs(scan.PrecursorCharge))
                {
                    if (includeExperimentalEnvelope)
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIonWithEnvelope(product, closestExperimentalMass.MonoisotopicMass.ToMz(closestExperimentalMass.Charge),
                            closestExperimentalMass.Peaks.First().intensity, closestExperimentalMass.Charge)
                        {
                            Envelope = closestExperimentalMass
                        });
                    }
                    else
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(product, closestExperimentalMass.MonoisotopicMass.ToMz(closestExperimentalMass.Charge),
                            closestExperimentalMass.Peaks.First().intensity, closestExperimentalMass.Charge));
                    }
                }
            }
            if (commonParameters.AddCompIons)
            {
                foreach (double massShift in ComplementaryIonConversion[commonParameters.DissociationType])
                {
                    double protonMassShift = massShift.ToMass(1);

                    for (int i = 0; i < theoreticalProducts.Count; i++)
                    {
                        var product = theoreticalProducts[i];
                        // unknown fragment mass or diagnostic ion or precursor; skip those
                        if (double.IsNaN(product.NeutralMass) || product.ProductType == ProductType.D || product.ProductType == ProductType.M)
                        {
                            continue;
                        }

                        double compIonMass = scan.PrecursorMass + protonMassShift - product.NeutralMass;

                        // get the closest peak in the spectrum to the theoretical peak
                        IsotopicEnvelope closestExperimentalMass = scan.GetClosestExperimentalIsotopicEnvelope(compIonMass);

                        // is the mass error acceptable?
                        if (closestExperimentalMass != null && productMassTolerance.Within(closestExperimentalMass.MonoisotopicMass, compIonMass) && closestExperimentalMass.Charge <= scan.PrecursorCharge)
                        {
                            //found the peak, but we don't want to save that m/z because it's the complementary of the observed ion that we "added". Need to create a fake ion instead.
                            double mz = (scan.PrecursorMass + protonMassShift - closestExperimentalMass.MonoisotopicMass).ToMz(closestExperimentalMass.Charge);

                            if (includeExperimentalEnvelope)
                            {
                                matchedFragmentIons.Add(new MatchedFragmentIonWithEnvelope(product, mz, closestExperimentalMass.TotalIntensity, closestExperimentalMass.Charge)
                                {
                                    Envelope = closestExperimentalMass
                                });
                            }
                            else
                            {
                                matchedFragmentIons.Add(new MatchedFragmentIon(product, mz, closestExperimentalMass.TotalIntensity, closestExperimentalMass.Charge));
                            }
                        }
                    }
                }
            }

            return matchedFragmentIons;
        }
    }
}
