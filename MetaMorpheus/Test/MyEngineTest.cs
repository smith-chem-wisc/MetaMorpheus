using EngineLayer;
using NUnit.Framework;
using System.Collections.Generic;
using System.Text;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;
using System;

namespace Test
{
    [TestFixture]
    public static class MyEngineTest
    {

        [Test]
        public static void TestMyEngine()
        {
            MetaMorpheusEngine level0engine = new TestEngine(0);

            level0engine = new TestEngine(0);
            level0engine.Run();
        }

        [Test]
        public static void MetaMorpheusEngineMatchFragmentIonsWithUnknownMass()
        {
            TestDataFile t = new TestDataFile();
            Tolerance productMassTolerance = new AbsoluteTolerance(0.01);
            double precursorMass = 300;
            //The below theoretical does not accurately represent B-Y ions
            double[] sorted_theoretical_product_masses_for_this_peptide = new double[] { precursorMass + (2 * Constants.ProtonMass) - 275.1350, precursorMass + (2 * Constants.ProtonMass) - 258.127, precursorMass + (2 * Constants.ProtonMass) - 257.1244, 50, 60, 70, 147.0764, precursorMass + (2 * Constants.ProtonMass) - 147.0764, precursorMass + (2 * Constants.ProtonMass) - 70, precursorMass + (2 * Constants.ProtonMass) - 60, precursorMass + (2 * Constants.ProtonMass) - 50, 257.1244, 258.127, 275.1350 }; //{ 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 }
            List<Product> productsWithLocalizedMassDiff = new List<Product>();
            foreach (double d in sorted_theoretical_product_masses_for_this_peptide)
            {
                productsWithLocalizedMassDiff.Add(new Product(ProductType.b, FragmentationTerminus.Both, d, 1, 1, 0));
            }

            Product productWithUnknownMass = new Product(ProductType.b, FragmentationTerminus.Both, Double.NaN, 1, 1, 0);
            productsWithLocalizedMassDiff.Add(productWithUnknownMass);

            CommonParameters commonParametersNoComp = new CommonParameters { ProductMassTolerance = new AbsoluteTolerance(0.01) };
            MsDataScan scan = t.GetOneBasedScan(2);

            //test xcorr processed spectrum with unknown mass (Double.NaN) which happens for unknown amino acid
            scan.MassSpectrum.XCorrPrePreprocessing(0, 1969, precursorMass.ToMz(1));
            var scanWithMass = new Ms2ScanWithSpecificMass(scan, precursorMass.ToMz(1), 1, "", new CommonParameters());
            List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, productsWithLocalizedMassDiff, commonParametersNoComp);
            Assert.That(1, Is.EqualTo(matchedIons.Count));

            //test what happens when the scan has no peaks
            scan.MassSpectrum.XCorrPrePreprocessing(0, 1, precursorMass.ToMz(1));
            scanWithMass = new Ms2ScanWithSpecificMass(scan, precursorMass.ToMz(1), 1, "", new CommonParameters());
            matchedIons = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, productsWithLocalizedMassDiff, commonParametersNoComp);
            Assert.That(0, Is.EqualTo(matchedIons.Count));
        }

[Test]
        public static void MatchFragmentIons_IncludeEnvelopeParameter()
        {
            TestDataFile t = new TestDataFile();
            double precursorMass = 402.18629720155; // Use the actual precursor mass from TestDataFile
            
            // Create theoretical products that match the spectrum peaks
            // TestDataFile scan2 has m/z peaks: { 50, 60, 70,147.0764, 257.1244, 258.127, 275.1350 }
            // The spectrum has isotopic envelopes when not XCorr processed
            List<Product> products = new List<Product>
            {
                // Monoisotopic neutral mass for m/z147.0764 at charge 1:147.0764 -1.007276 ≈ 146.069
                new Product(ProductType.b, FragmentationTerminus.Both,147.0764 - Constants.ProtonMass, 1, 1, 0)
            };

            CommonParameters commonParameters = new CommonParameters { ProductMassTolerance = new AbsoluteTolerance(0.5) };
            MsDataScan scan = t.GetOneBasedScan(2);
            
            // DO NOT call XCorrPrePreprocessing - we want to use the isotopic envelope path
            // which is only available when XcorrProcessed is false
            
            var scanWithMass = new Ms2ScanWithSpecificMass(scan, precursorMass.ToMz(2), 2, "", new CommonParameters());

            // Test default (no envelope)
            var matchedIonsDefault = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, products, commonParameters);
            
            // Test with envelope = false (explicit)
            var matchedIonsNoEnvelope = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, products, commonParameters, includeExperimentalEnvelope: false);
            
            // Test with envelope = true
            var matchedIonsWithEnvelope = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, products, commonParameters, includeExperimentalEnvelope: true);

            // Verify all return results
            Assert.That(matchedIonsDefault.Count, Is.GreaterThan(0), "Should match at least one ion");

            // Verify counts are same for all variants
            Assert.That(matchedIonsNoEnvelope.Count, Is.EqualTo(matchedIonsDefault.Count));
            Assert.That(matchedIonsWithEnvelope.Count, Is.EqualTo(matchedIonsDefault.Count));

            // Verify type when includeEnvelope = true
            Assert.That(matchedIonsWithEnvelope[0], Is.InstanceOf<MatchedFragmentIonWithEnvelope>());
            var ionWithEnvelope = (MatchedFragmentIonWithEnvelope)matchedIonsWithEnvelope[0];
            Assert.That(ionWithEnvelope.Envelope, Is.Not.Null);
            Assert.That(ionWithEnvelope.Envelope.Peaks, Is.Not.Empty);

            // Verify type when includeEnvelope = false
Assert.That(matchedIonsDefault[0], Is.InstanceOf<MatchedFragmentIon>());
            Assert.That(matchedIonsDefault[0], Is.Not.InstanceOf<MatchedFragmentIonWithEnvelope>());
        }

        [Test]
        public static void MatchFragmentIons_IncludeEnvelopeParameter_WithCompIons()
        {
            TestDataFile t = new TestDataFile();
            double precursorMass = 402.18629720155;
            
            List<Product> products = new List<Product>
            {
                new Product(ProductType.b, FragmentationTerminus.Both, 147.0764 - Constants.ProtonMass, 1, 1, 0)
            };

            CommonParameters commonParametersWithComp = new CommonParameters
            {
                ProductMassTolerance = new AbsoluteTolerance(0.5),
                AddCompIons = true
            };
            
            MsDataScan scan = t.GetOneBasedScan(2);
            var scanWithMass = new Ms2ScanWithSpecificMass(scan, precursorMass.ToMz(2), 2, "", new CommonParameters());

            // Test with complementary ions enabled and envelope = false
            var matchedIonsNoEnvelope = MetaMorpheusEngine.MatchFragmentIons(
                scanWithMass, products, commonParametersWithComp, matchAllCharges: false, includeExperimentalEnvelope: false);

            // Test with complementary ions enabled and envelope = true
            var matchedIonsWithEnvelope = MetaMorpheusEngine.MatchFragmentIons(
                scanWithMass, products, commonParametersWithComp, matchAllCharges: false, includeExperimentalEnvelope: true);

            // Both should find ions (main + complementary)
            Assert.That(matchedIonsNoEnvelope.Count, Is.GreaterThan(0), "Should match ions without envelope");
            Assert.That(matchedIonsWithEnvelope.Count, Is.EqualTo(matchedIonsNoEnvelope.Count), "Counts should match");

            // Verify envelope types - at least one should be a complementary ion if matched
            foreach (var ion in matchedIonsWithEnvelope)
            {
                if (ion is MatchedFragmentIonWithEnvelope ionWithEnv)
                {
                    Assert.That(ionWithEnv.Envelope, Is.Not.Null, "Envelope property should be populated");
                    Assert.That(ionWithEnv.Envelope.Peaks, Is.Not.Empty, "Envelope should have peaks");
                }
            }
        }

        private class TestEngine : MetaMorpheusEngine
        {

            public TestEngine(int level) : base(new CommonParameters(), null, new List<string>())
            {
            }

            protected override MetaMorpheusEngineResults RunSpecific()
            {
                return new TestResults(this);
            }

            private class TestResults : MetaMorpheusEngineResults
            {

                public TestResults(MetaMorpheusEngine e) : base(e)
                {
                }

                public override string ToString()
                {
                    var sb = new StringBuilder();
                    sb.AppendLine(base.ToString());
                    sb.Append("String for the TestResults results class");
                    return sb.ToString();
                }

            }

        }

    }
}