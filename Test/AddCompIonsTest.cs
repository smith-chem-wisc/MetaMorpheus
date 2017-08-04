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
    internal class AddCompIonsTest
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
            var searchModes = new List<MassDiffAcceptor> { new OpenSearchMode() };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            int maximumMissedCleavages = 0;
            int? minPeptideLength = null;
            int? maxPeptideLength = null;
            int maximumVariableModificationIsoforms = 4096;
            var engine = new ClassicSearchEngine(listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, minPeptideLength, maxPeptideLength, maximumVariableModificationIsoforms, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), false, InitiatorMethionineBehavior.Variable, false);
            var searchResults = (SearchResults)engine.Run();
            var engine2 = new ClassicSearchEngine(listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, minPeptideLength, maxPeptideLength, maximumVariableModificationIsoforms, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), false, InitiatorMethionineBehavior.Variable, true);
            var searchResults2 = (SearchResults)engine2.Run();

            // Single search mode
            Assert.AreEqual(searchResults.Psms.Length, searchResults2.Psms.Length);

            // Single ms2 scan
            Assert.AreEqual(searchResults.Psms[0].Length, searchResults2.Psms[0].Length);

            Assert.IsTrue(searchResults2.Psms[0][0].Score > 1);

            Assert.AreEqual(searchResults.Psms[0][0].ScanNumber, searchResults2.Psms[0][0].ScanNumber);

            Assert.IsTrue(searchResults2.Psms[0][0].Score < searchResults.Psms[0][0].Score * 2 && searchResults2.Psms[0][0].Score + 2 > searchResults.Psms[0][0].Score * 2);
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

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5) };
            var protease = new Protease("singleN", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, protease, initiatorMethionineBehavior, 2, null, null, 4096, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>());
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

            var engine = new ModernSearchEngine(listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, productMassTolerance, searchModes, new List<string>(), false, new List<ProductType> { ProductType.B, ProductType.Y });
            var searchResults = (SearchResults)engine.Run();
            var engine2 = new ModernSearchEngine(listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, productMassTolerance, searchModes, new List<string>(), true, new List<ProductType> { ProductType.B, ProductType.Y });
            var searchResults2 = (SearchResults)engine2.Run();

            // Single search mode
            Assert.AreEqual(searchResults.Psms.Length, searchResults2.Psms.Length);

            // Single ms2 scan
            Assert.AreEqual(searchResults.Psms[0].Length, searchResults2.Psms[0].Length);

            Assert.IsTrue(searchResults2.Psms[0][0].Score > 1);

            Assert.AreEqual(searchResults.Psms[0][0].ScanNumber, searchResults2.Psms[0][0].ScanNumber);

            Assert.IsTrue(searchResults2.Psms[0][0].Score < searchResults.Psms[0][0].Score * 2 && searchResults2.Psms[0][0].Score + 2 > searchResults.Psms[0][0].Score * 2);
        }

        [Test]
        public static void TestCompIons_MatchIons()
        {
            TestDataFile t = new TestDataFile();
            Tolerance productMassTolerance = new AbsoluteTolerance(0.01);
            double precursorMass = 300;
            double[] sorted_theoretical_product_masses_for_this_peptide = new double[] { precursorMass + (2 * Constants.protonMass) - 275.1350, precursorMass + (2 * Constants.protonMass) - 258.127, precursorMass + (2 * Constants.protonMass) - 257.1244, 50, 60, 70, 147.0764, precursorMass + (2 * Constants.protonMass) - 147.0764, precursorMass + (2 * Constants.protonMass) - 70, precursorMass + (2 * Constants.protonMass) - 60, precursorMass + (2 * Constants.protonMass) - 50, 257.1244, 258.127, 275.1350 }; //{ 50, 60, 70, 147.0764, 257.1244, 258.127, 275.1350 }
            double[] matchedIonMassesListPositiveIsMatch = new double[sorted_theoretical_product_masses_for_this_peptide.Count()];
            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            double scoreT = Psm.MatchIons(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, matchedIonMassesListPositiveIsMatch, true, precursorMass, lp);
            double scoreF = Psm.MatchIons(t.GetOneBasedScan(2), productMassTolerance, sorted_theoretical_product_masses_for_this_peptide, matchedIonMassesListPositiveIsMatch, false, precursorMass, lp);
            Assert.IsTrue(scoreT < scoreF * 2 && scoreT + 1 > scoreF * 2);
        }

        #endregion Public Methods

    }
}