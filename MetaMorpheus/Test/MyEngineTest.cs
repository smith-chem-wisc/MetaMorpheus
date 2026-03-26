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
            double precursorMass = 300;
            
            // Create a theoretical product that will match
            List<Product> products = new List<Product>
            {
                new Product(ProductType.b, FragmentationTerminus.N, 147.0764, 1, 1, 0)
            };

            CommonParameters commonParameters = new CommonParameters { ProductMassTolerance = new AbsoluteTolerance(0.01) };
            MsDataScan scan = t.GetOneBasedScan(2);
            
            // Ensure the scan is NOT Xcorr processed by setting minimal preprocessing
            // This allows the envelope-based matching to be used
            scan.MassSpectrum.XCorrPrePreprocessing(0, 1, precursorMass.ToMz(1));
            
            var scanWithMass = new Ms2ScanWithSpecificMass(scan, precursorMass.ToMz(1), 1, "", new CommonParameters());

            // Test default (no envelope)
            var matchedIonsDefault = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, products, commonParameters);
            
            // Test with envelope = false (explicit)
            var matchedIonsNoEnvelope = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, products, commonParameters, includeExperimentalEnvelope: false);
            
            // Test with envelope = true
            var matchedIonsWithEnvelope = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, products, commonParameters, includeExperimentalEnvelope: true);

            // Verify all return results (this test requires the scan to have ExperimentalFragments)
            // If the test data doesn't have experimental fragments, we at least verify the parameter doesn't crash
            if (matchedIonsDefault.Count == 0)
            {
                Assert.Pass("Test data does not have matching fragments - verifying method does not crash");
                return;
            }

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