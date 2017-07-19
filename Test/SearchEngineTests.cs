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
    public class SearchEngineTests
    {

        #region Public Methods

        [Test]
        public static void TestClassicSearchEngine()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5) };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            int maximumMissedCleavages = 2;
            int? minPeptideLength = null;
            int? maxPeptideLength = null;
            int maximumVariableModificationIsoforms = 4096;
            var engine = new ClassicSearchEngine(listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, minPeptideLength, maxPeptideLength, maximumVariableModificationIsoforms, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), false, InitiatorMethionineBehavior.Variable);
            var searchResults = (SearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.Psms.Length);

            // One scan
            Assert.AreEqual(1, searchResults.Psms[0].Length);

            Assert.IsTrue(searchResults.Psms[0][0].Score > 1);
            Assert.AreEqual(2, searchResults.Psms[0][0].ScanNumber);

            new SequencesToActualProteinPeptidesEngine(new List<SingleScanManyPeptidesMatch>[] { new List<SingleScanManyPeptidesMatch> { searchResults.Psms[0][0] } }, proteinList, searchModes, protease, maximumMissedCleavages, null, null, InitiatorMethionineBehavior.Variable, fixedModifications, variableModifications, 4096, new List<string>()).Run();

            Assert.AreEqual("QQQ", searchResults.Psms[0][0].MostProbableProteinInfo.BaseSequence);
        }

        [Test]
        public static void TestClassicSearchEngineWithWeirdPeptide()
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
            var engine = new ClassicSearchEngine(listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, minPeptideLength, maxPeptideLength, maximumVariableModificationIsoforms, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), false, InitiatorMethionineBehavior.Variable);
            var searchResults = (SearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.Psms.Length);

            // One Scan
            Assert.AreEqual(1, searchResults.Psms[0].Length);

            Assert.IsTrue(searchResults.Psms[0][0].Score > 1);
            Assert.AreEqual(2, searchResults.Psms[0][0].ScanNumber);

            new SequencesToActualProteinPeptidesEngine(new List<SingleScanManyPeptidesMatch>[] { new List<SingleScanManyPeptidesMatch> { searchResults.Psms[0][0] } }, proteinList, searchModes, protease, maximumMissedCleavages, null, null, InitiatorMethionineBehavior.Variable, fixedModifications, variableModifications, 4096, new List<string>()).Run();

            Assert.AreEqual("QXQ", searchResults.Psms[0][0].MostProbableProteinInfo.BaseSequence);
        }

        [Test]
        public static void TestModernSearchEngine()
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
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

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

            var engine = new ModernSearchEngine(listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, productMassTolerance, searchModes, new List<string>());
            var searchResults = (SearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.Psms.Length);

            // Single ms2 scan
            Assert.AreEqual(1, searchResults.Psms[0].Length);

            Assert.IsTrue(searchResults.Psms[0][0].Score > 1);
            Assert.AreEqual(2, searchResults.Psms[0][0].ScanNumber);

            new SequencesToActualProteinPeptidesEngine(new List<SingleScanManyPeptidesMatch>[] { new List<SingleScanManyPeptidesMatch> { searchResults.Psms[0][0] } }, proteinList, searchModes, protease, 2, null, null, initiatorMethionineBehavior, fixedModifications, variableModifications, 4096, new List<string>()).Run();

            Assert.AreEqual("QQQ", searchResults.Psms[0][0].MostProbableProteinInfo.BaseSequence);
        }

        [Test]
        public static void TestModernSearchEngineWithWeirdPeptide()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();

            var proteinList = new List<Protein> { new Protein("MNNNKQXQ", null) };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new List<MassDiffAcceptor> { new OpenSearchMode() };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            int maximumMissedCleavages = 2;
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, protease, initiatorMethionineBehavior, maximumMissedCleavages, null, null, 4096, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>());
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

            var engine = new ModernSearchEngine(listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, productMassTolerance, searchModes, new List<string>());
            var searchResults = (SearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.Psms.Length);

            // Single ms2 scan
            Assert.AreEqual(1, searchResults.Psms[0].Length);

            Assert.IsTrue(searchResults.Psms[0][0].Score > 1);
            Assert.AreEqual(2, searchResults.Psms[0][0].ScanNumber);

            new SequencesToActualProteinPeptidesEngine(new List<SingleScanManyPeptidesMatch>[] { new List<SingleScanManyPeptidesMatch> { searchResults.Psms[0][0] } }, proteinList, searchModes, protease, maximumMissedCleavages, null, null, initiatorMethionineBehavior, fixedModifications, variableModifications, 4096, new List<string>()).Run();

            Assert.AreEqual(3, searchResults.Psms[0][0].NumAmbiguous);
        }

        #endregion Public Methods

    }
}