using EngineLayer;
using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.Collections.Generic;
using System.Text;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;

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
        public static void MetaMorpheusEngineMatchFragmentIons()
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
            CommonParameters commonParametersNoComp = new CommonParameters { ProductMassTolerance = new AbsoluteTolerance(0.01) };
            CommonParameters commonParametersWithComp = new CommonParameters(productMassTolerance: new AbsoluteTolerance(0.01), addCompIons: true);

            MsDataScan scan = t.GetOneBasedScan(2);
            var scanWithMass = new Ms2ScanWithSpecificMass(scan, precursorMass.ToMz(1), 1, "", new CommonParameters());
            List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, productsWithLocalizedMassDiff, commonParametersNoComp);
             
            List<MatchedFragmentIon> matchedCompIons = MetaMorpheusEngine.MatchFragmentIons(scanWithMass, productsWithLocalizedMassDiff, commonParametersWithComp);
            matchedCompIons.AddRange(matchedIons);

            // score when the mass-diff is on this residue
            double localizedScore = MetaMorpheusEngine.CalculatePeptideScore(scan, matchedIons);
            double scoreNormal = MetaMorpheusEngine.CalculatePeptideScore(scan, matchedIons);
            double scoreComp = MetaMorpheusEngine.CalculatePeptideScore(scan, matchedCompIons);
            Assert.IsTrue(scoreNormal * 2 == scoreComp && scoreComp > scoreNormal + 1);
            
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