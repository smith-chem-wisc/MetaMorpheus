﻿using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
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
            Psm[][] allPsmsArray = new Psm[searchModes.Count()][];
            for (int aede = 0; aede < searchModes.Count; aede++)
                allPsmsArray[aede] = new Psm[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, minPeptideLength, maxPeptideLength, maximumVariableModificationIsoforms, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), false, InitiatorMethionineBehavior.Variable, false,1).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // One scan
            Assert.AreEqual(1, allPsmsArray[0].Length);

            Assert.IsTrue(allPsmsArray[0][0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0][0].ScanNumber);

            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm>[] { new List<Psm> { allPsmsArray[0][0] } }, proteinList, searchModes, protease, maximumMissedCleavages, null, null, InitiatorMethionineBehavior.Variable, fixedModifications, variableModifications, 4096, new List<string>(), TerminusType.None).Run();

            foreach (var huh in allPsmsArray[0])
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

            Assert.AreEqual("QQQ", allPsmsArray[0][0].BaseSequence);
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
            Psm[][] allPsmsArray = new Psm[searchModes.Count()][];
            for (int aede = 0; aede < searchModes.Count; aede++)
                allPsmsArray[aede] = new Psm[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, productMassTolerance, protease, searchModes, maximumMissedCleavages, minPeptideLength, maxPeptideLength, maximumVariableModificationIsoforms, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), false, InitiatorMethionineBehavior.Variable, false, 1).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // One Scan
            Assert.AreEqual(1, allPsmsArray[0].Length);

            Assert.IsTrue(allPsmsArray[0][0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0][0].ScanNumber);

            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm>[] { new List<Psm> { allPsmsArray[0][0] } }, proteinList, searchModes, protease, maximumMissedCleavages, null, null, InitiatorMethionineBehavior.Variable, fixedModifications, variableModifications, 4096, new List<string>(), TerminusType.None).Run();

            foreach (var huh in allPsmsArray[0])
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

            Assert.AreEqual("QXQ", allPsmsArray[0][0].BaseSequence);
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
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, protease, initiatorMethionineBehavior, 2, null, null, 4096, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), 1, 1, true);
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

            Psm[][] allPsmsArray = new Psm[searchModes.Count()][];
            for (int aede = 0; aede < searchModes.Count; aede++)
                allPsmsArray[aede] = new Psm[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, productMassTolerance, searchModes, new List<string>(), false, new List<ProductType>(), 1, 0, 1).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray[0].Length);

            Assert.IsTrue(allPsmsArray[0][0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0][0].ScanNumber);

            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm>[] { new List<Psm> { allPsmsArray[0][0] } }, proteinList, searchModes, protease, 2, null, null, initiatorMethionineBehavior, fixedModifications, variableModifications, 4096, new List<string>(), TerminusType.None).Run();

            foreach (var huh in allPsmsArray[0])
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

            Assert.AreEqual("QQQ", allPsmsArray[0][0].BaseSequence);
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
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, protease, initiatorMethionineBehavior, maximumMissedCleavages, null, null, 4096, new List<ProductType> { ProductType.B, ProductType.Y }, new List<string>(), 1, 1, true);
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

            Psm[][] allPsmsArray = new Psm[searchModes.Count()][];
            for (int aede = 0; aede < searchModes.Count; aede++)
                allPsmsArray[aede] = new Psm[listOfSortedms2Scans.Length];
            var engine = new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, productMassTolerance, searchModes, new List<string>(), false, new List<ProductType>(), 1, 0, 1);
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray[0].Length);

            Assert.IsTrue(allPsmsArray[0][0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0][0].ScanNumber);

            new SequencesToActualProteinPeptidesEngine(new List<Psm>[] { new List<Psm> { allPsmsArray[0][0] } }, proteinList, searchModes, protease, maximumMissedCleavages, null, null, initiatorMethionineBehavior, fixedModifications, variableModifications, 4096, new List<string>(), TerminusType.None).Run();

            Assert.AreEqual(3, allPsmsArray[0][0].NumDifferentCompactPeptides);
        }

        [Test]
        public static void TestNonSpecificEnzymeEngineSingleN()
        {
            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
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

            var proteinList = new List<Protein> { new Protein("GGGGGMNNNKQQQGGGGG", null) };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5), new OpenSearchMode() };
            var protease = new Protease("singleN", new List<string> { "K, G" }, new List<string>(), TerminusType.None, CleavageSpecificity.None, null, null, null);

            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, protease, initiatorMethionineBehavior, 2, null, null, 4096, new List<ProductType> { ProductType.B }, new List<string>(), 1, 1, true);
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

            Psm[][] allPsmsArray = new Psm[searchModes.Count()][];
            for (int aede = 0; aede < searchModes.Count; aede++)
                allPsmsArray[aede] = new Psm[listOfSortedms2Scans.Length];
            var engine = new NonSpecificEnzymeEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, productMassTolerance, searchModes, new List<string>(), true, new List<ProductType> { ProductType.B }, protease, 5, TerminusType.N, 1, 0, 1);
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(2, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray[1].Length);

            Assert.IsTrue(allPsmsArray[1][0].Score > 4);
            Assert.AreEqual(2, allPsmsArray[1][0].ScanNumber);

            var hah = (SequencesToActualProteinPeptidesEngineResults)new NonSpecificEnzymeSequencesToActualPeptides(new List<Psm>[] { new List<Psm> { allPsmsArray[1][0] } }, proteinList, searchModes, protease, 2, null, null, initiatorMethionineBehavior, fixedModifications, variableModifications, 4096, new List<string>(), TerminusType.N).Run();

            foreach (var huh in allPsmsArray[1])
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

            Assert.AreEqual("QQQGGGG", allPsmsArray[1][0].BaseSequence);
        }

        [Test]
        public static void TestNonSpecificEnzymeEngineSingleC()
        {
            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
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

            var proteinList = new List<Protein> { new Protein("GGGGGMNNNKQQQGGGGG", null) };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5), new OpenSearchMode() };
            var protease = new Protease("singleC", new List<string> { "K, G" }, new List<string>(), TerminusType.None, CleavageSpecificity.None, null, null, null);

            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, protease, initiatorMethionineBehavior, 2, null, null, 4096, new List<ProductType> { ProductType.Y }, new List<string>(), 1, 1, true);
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

            Psm[][] allPsmsArray = new Psm[searchModes.Count()][];
            for (int aede = 0; aede < searchModes.Count; aede++)
                allPsmsArray[aede] = new Psm[listOfSortedms2Scans.Length];
            var engine = new NonSpecificEnzymeEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, productMassTolerance, searchModes, new List<string>(), true, new List<ProductType> { ProductType.Y }, protease, 5, TerminusType.C, 1, 0, 1);
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(2, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray[1].Length);

            Assert.IsTrue(allPsmsArray[1][0].Score > 4);
            Assert.AreEqual(2, allPsmsArray[1][0].ScanNumber);

            var hah = (SequencesToActualProteinPeptidesEngineResults)new NonSpecificEnzymeSequencesToActualPeptides(new List<Psm>[] { new List<Psm> { allPsmsArray[1][0] } }, proteinList, searchModes, protease, 2, null, null, initiatorMethionineBehavior, fixedModifications, variableModifications, 4096, new List<string>(), TerminusType.C).Run();

            foreach (var huh in allPsmsArray[1])
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

            Assert.AreEqual("QQQGGGG", allPsmsArray[1][0].BaseSequence);
        }

        #endregion Public Methods
    }
}