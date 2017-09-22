using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    internal static class AddCompIonsTest
    {
        #region Public Methods

        [Test]
        public static void TestAddCompIonsClassic()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("QXQ", null) };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new OpenSearchMode();

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];

            CommonParameters CommonParameters = new CommonParameters
            {
                DigestionParams = new DigestionParams
                {
                    Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
                    MinPeptideLength = null,
                    MaxMissedCleavages = 0
                },
                ConserveMemory = false,
                ScoreCutoff = 1,
            };
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, false, CommonParameters, new List<string>()).Run();
            Psm[] allPsmsArray2 = new Psm[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray2, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, true, CommonParameters, new List<string>()).Run();

            double scoreT = allPsmsArray2[0].Score;
            double scoreF = allPsmsArray[0].Score;

            // Single search mode
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            // Single ms2 scan
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            Assert.IsTrue(scoreT > 1);

            Assert.AreEqual(allPsmsArray[0].ScanNumber, allPsmsArray2[0].ScanNumber);

            Assert.IsTrue(scoreT == scoreF * 3 && scoreT > scoreF + 2);
        }

        [Test]
        public static void TestCompIons_ModernSearch()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int ii = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)ii);
                ii++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)ii);
                ii++;
            }

            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            SearchParameters SearchParameters = new SearchParameters
            {
                MassDiffAcceptor = searchModes
            };
            CommonParameters CommonParameters = new CommonParameters
            {
                DigestionParams = new DigestionParams
                {
                    Protease = new Protease("singleN", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
                    MinPeptideLength = null,
                },
                ConserveMemory = false,
                ScoreCutoff = 1,
            };
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B, ProductType.Y }, 1, true, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.TotalPartitions, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, SearchParameters.AddCompIons, SearchParameters.MassDiffAcceptor, new List<string>()).Run();

            Psm[] allPsmsArray2 = new Psm[listOfSortedms2Scans.Length];
            SearchParameters.AddCompIons = true;
            new ModernSearchEngine(allPsmsArray2, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, SearchParameters.AddCompIons, SearchParameters.MassDiffAcceptor, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            // Single ms2 scan
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            Assert.IsTrue(allPsmsArray2[0].Score > 1);

            Assert.AreEqual(allPsmsArray[0].ScanNumber, allPsmsArray2[0].ScanNumber);

            Assert.IsTrue(allPsmsArray2[0].Score <= allPsmsArray[0].Score * 2 && allPsmsArray2[0].Score > allPsmsArray[0].Score + 3);
        }

        [Test]
        public static void TestCompIons_MatchIonsScore()
        {
            TestDataFile t = new TestDataFile();
            Tolerance productMassTolerance = new AbsoluteTolerance(0.01);
            double precursorMass = 300;
            //The below theoretical does not accurately represent B-Y ions
            double[] sorted_theoretical_product_masses_for_this_peptide = new double[] { precursorMass + (2 * Constants.protonMass) - 275.1350, precursorMass + (2 * Constants.protonMass) - 258.127, precursorMass + (2 * Constants.protonMass) - 257.1244, 50, 60, 70, 147.0764, precursorMass + (2 * Constants.protonMass) - 147.0764, precursorMass + (2 * Constants.protonMass) - 70, precursorMass + (2 * Constants.protonMass) - 60, precursorMass + (2 * Constants.protonMass) - 50, 257.1244, 258.127, 275.1350 }; //{ 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 }
            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            double scoreT = MetaMorpheusEngine.CalculatePeptideScore(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, precursorMass, new List<DissociationType> { DissociationType.HCD }, true);
            double scoreF = MetaMorpheusEngine.CalculatePeptideScore(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, precursorMass, new List<DissociationType> { DissociationType.HCD }, false);
            Assert.IsTrue(scoreT == scoreF * 2 && scoreT > scoreF + 1);
        }

        [Test]
        public static void TestCompIons_MatchIons()
        {
            TestDataFile t = new TestDataFile(0.0001);
            Tolerance productMassTolerance = new AbsoluteTolerance(0.01);
            double precursorMass = 402.18629720155;
            //The below theoretical does not accurately represent B-Y ions
            double[] sorted_theoretical_product_masses_for_this_peptide = new double[] { 50, 60, 70, 147.0764 - Constants.protonMass, 200, 215, 230, 245, precursorMass + Constants.protonMass - 147.0764, 258.127, 275.1350, precursorMass + (2 * Constants.protonMass) - 70, precursorMass + (2 * Constants.protonMass) - 60, precursorMass + (2 * Constants.protonMass) - 50 }; //{ 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 }
            List<double> matchedIonsT = new List<double>();
            List<double> matchedDaErrorT = new List<double>();
            List<double> matchedPpmErrorT = new List<double>();
            List<double> matchedIonsF = new List<double>();
            List<double> matchedDaErrorF = new List<double>();
            List<double> matchedPpmErrorF = new List<double>();
            MetaMorpheusEngine.MatchIons(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, matchedIonsT, matchedDaErrorT, matchedPpmErrorT, precursorMass, new List<DissociationType> { DissociationType.HCD }, true);
            MetaMorpheusEngine.MatchIons(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, matchedIonsF, matchedDaErrorF, matchedPpmErrorF, precursorMass, new List<DissociationType> { DissociationType.HCD }, false);

            //Test the number of ions is doubled
            Assert.IsTrue(matchedIonsT.Count == matchedIonsF.Count * 2);
            //Test the number of da errors is doubled
            Assert.IsTrue(matchedDaErrorT.Count == matchedDaErrorF.Count * 2);
            //test the number of ppm errors is doubled
            Assert.IsTrue(matchedPpmErrorT.Count == matchedPpmErrorF.Count * 2);
            foreach (double d in matchedDaErrorF)
                Assert.IsTrue(d <= 0.01);
            foreach (double d in matchedDaErrorT)
                Assert.IsTrue(d <= 0.01);
        }

        #endregion Public Methods
    }
}