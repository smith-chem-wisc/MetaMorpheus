using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
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

            // Single search mode
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            // Single ms2 scan
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            Assert.IsTrue(allPsmsArray2[0].Score > 1);

            Assert.AreEqual(allPsmsArray[0].ScanNumber, allPsmsArray2[0].ScanNumber);

            Assert.IsTrue(allPsmsArray2[0].Score < allPsmsArray[0].Score * 2 && allPsmsArray2[0].Score + 2 > allPsmsArray[0].Score * 2);
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
                MassDiffAcceptors = searchModes
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
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndexDict;
            var keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
            var fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, SearchParameters.AddCompIons, SearchParameters.MassDiffAcceptors, new List<string>()).Run();

            Psm[] allPsmsArray2 = new Psm[listOfSortedms2Scans.Length];
            SearchParameters.AddCompIons = true;
            new ModernSearchEngine(allPsmsArray2, listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, SearchParameters.AddCompIons, SearchParameters.MassDiffAcceptors, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            // Single ms2 scan
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            Assert.IsTrue(allPsmsArray2[0].Score > 1);

            Assert.AreEqual(allPsmsArray[0].ScanNumber, allPsmsArray2[0].ScanNumber);

            Assert.IsTrue(allPsmsArray2[0].Score < allPsmsArray[0].Score * 2 && allPsmsArray2[0].Score + 2 > allPsmsArray[0].Score * 2);
        }

        [Test]
        public static void TestCompIons_MatchIons()
        {
            TestDataFile t = new TestDataFile();
            Tolerance productMassTolerance = new AbsoluteTolerance(0.01);
            double precursorMass = 300;
            double[] sorted_theoretical_product_masses_for_this_peptide = new double[] { precursorMass + (2 * Constants.protonMass) - 275.1350, precursorMass + (2 * Constants.protonMass) - 258.127, precursorMass + (2 * Constants.protonMass) - 257.1244, 50, 60, 70, 147.0764, precursorMass + (2 * Constants.protonMass) - 147.0764, precursorMass + (2 * Constants.protonMass) - 70, precursorMass + (2 * Constants.protonMass) - 60, precursorMass + (2 * Constants.protonMass) - 50, 257.1244, 258.127, 275.1350 }; //{ 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 }
            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            double scoreT = MetaMorpheusEngine.CalculateClassicScore(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, true, precursorMass, lp);
            double scoreF = MetaMorpheusEngine.CalculateClassicScore(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, false, precursorMass, lp);
            Assert.IsTrue(scoreT < scoreF * 2 && scoreT + 1 > scoreF * 2);
        }

        #endregion Public Methods
    }
}