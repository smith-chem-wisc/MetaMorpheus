﻿using EngineLayer;
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
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null, null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], null, null, false, false, null) };

            var productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            var searchModes = new List<SearchMode> { new SinglePpmAroundZeroSearchMode(5) };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var listOfSortedms2Scans = MetaMorpheusEngine.GetMs2Scans(myMsDataFile).OrderBy(b => b.MonoisotopicPrecursorMass).ToArray();
            int maximumMissedCleavages = 2;
            int? minPeptideLength = null;
            int? maxPeptideLength = null;
            int maximumVariableModificationIsoforms = 4096;
            var engine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, minPeptideLength, maxPeptideLength, maximumVariableModificationIsoforms, "lawl", new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), false);
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
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("MNNNKQXQ", null, null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], null, null, false, false, null) };

            var productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            var searchModes = new List<SearchMode> { new OpenSearchMode() };
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var listOfSortedms2Scans = MetaMorpheusEngine.GetMs2Scans(myMsDataFile).OrderBy(b => b.MonoisotopicPrecursorMass).ToArray();

            int maximumMissedCleavages = 2;
            int? minPeptideLength = null;
            int? maxPeptideLength = null;
            int maximumVariableModificationIsoforms = 4096;
            var engine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, minPeptideLength, maxPeptideLength, maximumVariableModificationIsoforms, "UUU", new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), false);
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

            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null, null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], null, null, false, false, null) };

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

            var engine = new ModernSearchEngine(myMsDataFile, peptideIndex, keys, fragmentIndex, productMassTolerance, searchModes, "fayk", null);
            var searchResults = (ModernSearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.NewPsms.Length);

            // Two scans, even including the MS1 scans
            Assert.AreEqual(2, searchResults.NewPsms[0].Count);

            Assert.IsTrue(searchResults.NewPsms[0][1].score > 1);
            Assert.AreEqual(2, searchResults.NewPsms[0][1].scanNumber);
            Assert.AreEqual("QQQ", searchResults.NewPsms[0][1].GetCompactPeptide(modsDictionary).BaseSequence);
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

            var proteinList = new List<Protein> { new Protein("MNNNKQXQ", null, null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], null, null, false, false, null) };

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

            var engine = new ModernSearchEngine(myMsDataFile, peptideIndex, keys, fragmentIndex, productMassTolerance, searchModes, "fayk", null);
            var searchResults = (ModernSearchResults)engine.Run();

            // Single search mode
            Assert.AreEqual(1, searchResults.NewPsms.Length);

            // Two scans, even including the MS1 scans
            Assert.AreEqual(2, searchResults.NewPsms[0].Count);

            Assert.IsTrue(searchResults.NewPsms[0][1].score > 1);
            Assert.AreEqual(2, searchResults.NewPsms[0][1].scanNumber);
            Assert.AreEqual("QXQ", searchResults.NewPsms[0][1].GetCompactPeptide(modsDictionary).BaseSequence);
        }

        #endregion Public Methods

    }
}