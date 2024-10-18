using EngineLayer;
using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
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
            Assert.AreEqual(1, matchedIons.Count);

            //test what happens when the scan has no peaks
            scan.MassSpectrum.XCorrPrePreprocessing(0, 1, precursorMass.ToMz(1));
            scanWithMass = new Ms2ScanWithSpecificMass(scan, precursorMass.ToMz(1), 1, "", new CommonParameters());
            matchedIons = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, productsWithLocalizedMassDiff, commonParametersNoComp);
            Assert.AreEqual(0, matchedIons.Count);
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