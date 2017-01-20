using InternalLogicEngineLayer;
using NUnit.Framework;
using OldInternalLogic;
using Spectra;
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
            var variableModifications = new List<MorpheusModification>();
            var fixedModifications = new List<MorpheusModification>();
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false) };

            var productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            var searchModes = new List<SearchMode> { new SinglePpmAroundZeroSearchMode("", 5) };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);

            var listOfSortedms2Scans = myMsDataFile.Where(b => b.MsnOrder == 2).Select(b => new LocalMS2Scan(b)).OrderBy(b => b.PrecursorMass).ToArray();

            int maximumMissedCleavages = 2;
            int maximumVariableModificationIsoforms = 4096;
            var engine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, maximumVariableModificationIsoforms, "lawl");
            var searchResults = (ClassicSearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.OuterPsms.Length);

            // Two scans, even including the MS1 scans
            Assert.AreEqual(2, searchResults.OuterPsms[0].Length);

            Assert.IsTrue(searchResults.OuterPsms[0][1].score > 1);
            Assert.AreEqual(2, searchResults.OuterPsms[0][1].scanNumber);
            Assert.AreEqual("QQQ", searchResults.OuterPsms[0][1].ps.BaseSequence);
        }

        [Test]
        public static void TestClassicSearchEngineWithWeirdPeptide()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<MorpheusModification>();
            var fixedModifications = new List<MorpheusModification>();
            var proteinList = new List<Protein> { new Protein("MNNNKQXQ", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false) };

            var productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            var searchModes = new List<SearchMode> { new OpenSearchMode("open") };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);

            var listOfSortedms2Scans = myMsDataFile.Where(b => b.MsnOrder == 2).Select(b => new LocalMS2Scan(b)).OrderBy(b => b.PrecursorMass).ToArray();

            int maximumMissedCleavages = 2;
            int maximumVariableModificationIsoforms = 4096;
            var engine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, maximumVariableModificationIsoforms, "UUU");
            var searchResults = (ClassicSearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.OuterPsms.Length);

            // Two scans, even including the MS1 scans
            Assert.AreEqual(2, searchResults.OuterPsms[0].Length);

            Assert.IsTrue(searchResults.OuterPsms[0][1].score > 1);
            Assert.AreEqual(2, searchResults.OuterPsms[0][1].scanNumber);
            Assert.AreEqual("QXQ", searchResults.OuterPsms[0][1].ps.BaseSequence);
        }

        [Test]
        public static void TestModernSearchEngine()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<MorpheusModification>();
            var fixedModifications = new List<MorpheusModification>();
            var localizeableModifications = new List<MorpheusModification>();
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false) };

            var productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            var searchModes = new List<SearchMode> { new SinglePpmAroundZeroSearchMode("", 5) };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);

            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            var indexEngine = new IndexEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, protease, initiatorMethionineBehavior, 2, 4096);
            var indexResults = (IndexResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndexDict;
            var keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
            var fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();

            var engine = new ModernSearchEngine(myMsDataFile, peptideIndex, keys, fragmentIndex, productMassTolerance.Value, searchModes);
            var searchResults = (ModernSearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.NewPsms.Length);

            // Two scans, even including the MS1 scans
            Assert.AreEqual(2, searchResults.NewPsms[0].Count);

            Assert.IsTrue(searchResults.NewPsms[0][1].score > 1);
            Assert.AreEqual(2, searchResults.NewPsms[0][1].scanNumber);
            Assert.AreEqual("QQQ", searchResults.NewPsms[0][1].GetCompactPeptide(variableModifications, localizeableModifications).BaseSequence);
        }

        [Test]
        public static void TestModernSearchEngineWithWeirdPeptide()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<MorpheusModification>();
            var fixedModifications = new List<MorpheusModification>();
            var localizeableModifications = new List<MorpheusModification>();
            var proteinList = new List<Protein> { new Protein("MNNNKQXQ", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false) };

            var productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            var searchModes = new List<SearchMode> { new OpenSearchMode("d") };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);

            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            var indexEngine = new IndexEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, protease, initiatorMethionineBehavior, 2, 4096);
            var indexResults = (IndexResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndexDict;
            var keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
            var fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();

            var engine = new ModernSearchEngine(myMsDataFile, peptideIndex, keys, fragmentIndex, productMassTolerance.Value, searchModes);
            var searchResults = (ModernSearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.NewPsms.Length);

            // Two scans, even including the MS1 scans
            Assert.AreEqual(2, searchResults.NewPsms[0].Count);

            Assert.IsTrue(searchResults.NewPsms[0][1].score > 1);
            Assert.AreEqual(2, searchResults.NewPsms[0][1].scanNumber);
            Assert.AreEqual("QXQ", searchResults.NewPsms[0][1].GetCompactPeptide(variableModifications, localizeableModifications).BaseSequence);
        }

        #endregion Public Methods

    }
}