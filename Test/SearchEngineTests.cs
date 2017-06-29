using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;

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

            var productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            var searchModes = new List<SearchMode> { new SinglePpmAroundZeroSearchMode(5) };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            bool useProvidedPrecursorInfo = true;
            bool findAllPrecursors = true;
            var intensityRatio = 4;
            var listOfSortedms2Scans = MetaMorpheusEngine.GetMs2Scans(myMsDataFile, findAllPrecursors, useProvidedPrecursorInfo, intensityRatio, null).OrderBy(b => b.PrecursorMass).ToArray();
            int maximumMissedCleavages = 2;
            int? minPeptideLength = null;
            int? maxPeptideLength = null;
            int maximumVariableModificationIsoforms = 4096;
            var engine = new ClassicSearchEngine(listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, minPeptideLength, maxPeptideLength, maximumVariableModificationIsoforms, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), false);
            var searchResults = (ClassicSearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.OuterPsms.Length);

            // One scan
            Assert.AreEqual(1, searchResults.OuterPsms[0].Length);

            Assert.IsTrue(searchResults.OuterPsms[0][0].Score > 1);
            Assert.AreEqual(2, searchResults.OuterPsms[0][0].ScanNumber);
            Assert.AreEqual("QQQ", (searchResults.OuterPsms[0][0] as PsmClassic).ps.BaseSequence);
        }

        [Test]
        public static void TestClassicSearchEngineWithWeirdPeptide()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("QXQ", null) };

            var productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            var searchModes = new List<SearchMode> { new OpenSearchMode() };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            bool useProvidedPrecursorInfo = true;
            bool findAllPrecursors = true;
            var intensityRatio = 4;
            var listOfSortedms2Scans = MetaMorpheusEngine.GetMs2Scans(myMsDataFile, findAllPrecursors, useProvidedPrecursorInfo, intensityRatio, null).OrderBy(b => b.PrecursorMass).ToArray();

            int maximumMissedCleavages = 0;
            int? minPeptideLength = null;
            int? maxPeptideLength = null;
            int maximumVariableModificationIsoforms = 4096;
            var engine = new ClassicSearchEngine(listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, minPeptideLength, maxPeptideLength, maximumVariableModificationIsoforms, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), false);
            var searchResults = (ClassicSearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.OuterPsms.Length);

            // One Scan
            Assert.AreEqual(1, searchResults.OuterPsms[0].Length);

            Assert.IsTrue(searchResults.OuterPsms[0][0].Score > 1);
            Assert.AreEqual(2, searchResults.OuterPsms[0][0].ScanNumber);
            Assert.AreEqual("QXQ", (searchResults.OuterPsms[0][0] as PsmClassic).ps.BaseSequence);
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

            var productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            var searchModes = new List<SearchMode> { new SinglePpmAroundZeroSearchMode(5) };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, modsDictionary, protease, initiatorMethionineBehavior, 2, null, null, 4096, new List<ProductType> { ProductType.B, ProductType.Y }, null);
            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndexDict;
            var keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
            var fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();
            bool useProvidedPrecursorInfo = true;
            bool findAllPrecursors = true;
            var intensityRatio = 4;
            var listOfSortedms2Scans = MetaMorpheusEngine.GetMs2Scans(myMsDataFile, findAllPrecursors, useProvidedPrecursorInfo, intensityRatio, null).OrderBy(b => b.PrecursorMass).ToArray();

            var engine = new ModernSearchEngine(listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, productMassTolerance, searchModes, null);
            var searchResults = (ModernSearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.OuterPsms.Length);

            // Single ms2 scan
            Assert.AreEqual(1, searchResults.OuterPsms[0].Length);

            Assert.IsTrue(searchResults.OuterPsms[0][0].Score > 1);
            Assert.AreEqual(2, searchResults.OuterPsms[0][0].ScanNumber);
            Assert.AreEqual("QQQ", searchResults.OuterPsms[0][0].GetCompactPeptide(modsDictionary).BaseSequence);
        }

        [Test]
        public static void TestModernSearchEngineWithWeirdPeptide()
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

            var proteinList = new List<Protein> { new Protein("MNNNKQXQ", null) };

            var productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            var searchModes = new List<SearchMode> { new OpenSearchMode() };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, modsDictionary, protease, initiatorMethionineBehavior, 2, null, null, 4096, new List<ProductType> { ProductType.B, ProductType.Y }, null);
            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndexDict;
            var keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
            var fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();
            bool useProvidedPrecursorInfo = true;
            bool findAllPrecursors = true;
            var intensityRatio = 4;
            var listOfSortedms2Scans = MetaMorpheusEngine.GetMs2Scans(myMsDataFile, findAllPrecursors, useProvidedPrecursorInfo, intensityRatio, null).OrderBy(b => b.PrecursorMass).ToArray();

            var engine = new ModernSearchEngine(listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, productMassTolerance, searchModes, null);
            var searchResults = (ModernSearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.OuterPsms.Length);

            // Single ms2 scan
            Assert.AreEqual(1, searchResults.OuterPsms[0].Length);

            Assert.IsTrue(searchResults.OuterPsms[0][0].Score > 1);
            Assert.AreEqual(2, searchResults.OuterPsms[0][0].ScanNumber);
            Assert.AreEqual("QXQ", searchResults.OuterPsms[0][0].GetCompactPeptide(modsDictionary).BaseSequence);
        }

        #endregion Public Methods

    }
}