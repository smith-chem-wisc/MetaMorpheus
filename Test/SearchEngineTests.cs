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
        public void TestClassicSearchEngine()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<MorpheusModification>();
            var fixedModifications = new List<MorpheusModification>();
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null, null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false) };

            var productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            var searchModes = new List<SearchMode> { new SinglePpmAroundZeroSearchMode("", 5) };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), Terminus.C, CleavageSpecificity.Full, null, null, null);

            var listOfSortedms2Scans = myMsDataFile.Where(b => b.MsnOrder == 2).Select(b => new LocalMs2Scan(b)).OrderBy(b => b.precursorMass).ToArray();

            var engine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, 0, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes);
            var searchResults = (ClassicSearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.outerPsms.Length);

            // Two scans, even including the MS1 scans
            Assert.AreEqual(2, searchResults.outerPsms[0].Length);

            Assert.IsTrue(searchResults.outerPsms[0][1].Score > 1);
            Assert.AreEqual(2, searchResults.outerPsms[0][1].scanNumber);
            Assert.AreEqual("QQQ", searchResults.outerPsms[0][1].ps.BaseSequence);
        }

        [Test]
        public void TestModernSearchEngine()
        {
            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<MorpheusModification>();
            var fixedModifications = new List<MorpheusModification>();
            var localizeableModifications = new List<MorpheusModification>();
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null, null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false) };

            var productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            var searchModes = new List<SearchMode> { new SinglePpmAroundZeroSearchMode("", 5) };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), Terminus.C, CleavageSpecificity.Full, null, null, null);

            var indexEngine = new IndexEngine(proteinList, variableModifications, fixedModifications, localizeableModifications, protease);
            var indexResults = (IndexResults)indexEngine.Run();
            var peptideIndex = indexResults.peptideIndex;
            var fragmentIndexDict = indexResults.fragmentIndexDict;
            var keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
            var fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();

            var engine = new ModernSearchEngine(myMsDataFile, 0, peptideIndex, keys, fragmentIndex, productMassTolerance.Value, searchModes);
            var searchResults = (ModernSearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.newPsms.Length);

            // Two scans, even including the MS1 scans
            Assert.AreEqual(2, searchResults.newPsms[0].Count);

            Assert.IsTrue(searchResults.newPsms[0][1].Score > 1);
            Assert.AreEqual(2, searchResults.newPsms[0][1].scanNumber);
            Assert.AreEqual("QQQ", searchResults.newPsms[0][1].GetCompactPeptide(variableModifications, localizeableModifications).BaseSequence);
        }

        #endregion Public Methods

    }
}