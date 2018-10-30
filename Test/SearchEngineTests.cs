﻿using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class SearchEngineTests
    {
        [Test]
        public static void TestClassicSearchEngine()
        {
            Protease protease = new Protease("Customized Protease", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters
                (digestionParams: new DigestionParams(
                    protease: protease.Name,
                    minPeptideLength: 1),
                scoreCutoff: 1);

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // One scan
            Assert.AreEqual(1, allPsmsArray.Length);

            Assert.IsTrue(allPsmsArray[0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

            Assert.AreEqual("QQQ", allPsmsArray[0].BaseSequence);
        }

        [Test]
        public static void TestClassicSearchEngineWithWeirdPeptide()
        {
            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(
                    protease: "Customized Protease",
                    maxMissedCleavages: 0,
                    minPeptideLength: 1),
                scoreCutoff: 1);

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("QXQ", null) };

            var searchModes = new OpenSearchMode();

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // One Scan
            Assert.AreEqual(1, allPsmsArray.Length);

            Assert.IsTrue(allPsmsArray[0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

            Assert.AreEqual("QXQ", allPsmsArray[0].BaseSequence);
        }

        [Test]
        public static void TestModernSearchEngine()
        {
            SearchParameters SearchParameters = new SearchParameters
            {
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                SearchTarget = true,
            };

            CommonParameters CommonParameters = new CommonParameters(
                precursorMassTolerance: new PpmTolerance(5),
                digestionParams: new DigestionParams(
                    protease: "Customized Protease",
                    minPeptideLength: 1),
                scoreCutoff: 1);

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();

            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, CommonParameters, SearchParameters.MaxFragmentSize, false, new List<FileInfo>(), new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray.Length);
            Assert.That(allPsmsArray[0] != null);

            Assert.IsTrue(allPsmsArray[0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

            Assert.AreEqual("QQQ", allPsmsArray[0].BaseSequence);
        }

        [Test]
        public static void TestModernSearchFragmentEdges()
        {
            // purpose of this test is to check that fragments at the edge of the fragment index in modern search do not fall out of bounds
            // example: fragment mass 400 when max frag mass is 400 may go over the edge of the array because of ppm tolerances
            SearchParameters SearchParameters = new SearchParameters
            {
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                MaxFragmentSize = 1, // super small index
            };
            CommonParameters CommonParameters = new CommonParameters(
                productMassTolerance: new AbsoluteTolerance(100), // super large tolerance (100 Da)
                digestionParams: new DigestionParams(
                    protease: "trypsin",
                    minPeptideLength: 1),
                scoreCutoff: 1,
                addCompIons: true);

            var proteinList = new List<Protein> { new Protein("K", null) };

            var indexEngine = new IndexingEngine(proteinList, new List<Modification>(), new List<Modification>(), 1, DecoyType.Reverse,
                CommonParameters, SearchParameters.MaxFragmentSize, false, new List<FileInfo>(), new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(new TestDataFile(), null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

            // no assertions... just don't crash...
        }

        [Test]
        public static void TestNonViablePSM()
        {
            //just check if it crashes or not.
            SearchParameters SearchParameters = new SearchParameters();

            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(
                    protease: "Customized Protease"),
                scoreCutoff: 255);

            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            var myMsDataFile = new TestDataFile(true); //empty
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();

            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, CommonParameters, SearchParameters.MaxFragmentSize, false, new List<FileInfo>(), new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[][] allPsmsArrays = new PeptideSpectralMatch[1][];
            allPsmsArrays[0] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            PeptideSpectralMatch[] allPsmsArray = allPsmsArrays[0];

            //Classic
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters, new List<string>()).Run();

            //Modern
            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

            //NonSpecific
            new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, indexResults.FragmentIndex, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();
        }

        [Test]
        public static void TestModernSearchEngineWithWeirdPeptide()
        {
            SearchParameters SearchParameters = new SearchParameters
            {
                MassDiffAcceptorType = MassDiffAcceptorType.Open,
                SearchTarget = true,
            };

            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(
                    protease: "Customized Protease",
                    minPeptideLength: 1),
                scoreCutoff: 1);

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();

            var proteinList = new List<Protein> { new Protein("MNNNKQXQ", null) };

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, CommonParameters, SearchParameters.MaxFragmentSize, false, new List<FileInfo>(), new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            var engine = new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>());
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray.Length);
            Assert.That(allPsmsArray[0] != null);

            Assert.IsTrue(allPsmsArray[0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

            Assert.AreEqual(3, allPsmsArray[0].NumDifferentMatchingPeptides);
        }

        [Test]
        public static void TestNonSpecificEnzymeSearchEngineSingleN()
        {
            SearchParameters SearchParameters = new SearchParameters
            {
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                LocalFdrCategories = new List<FdrCategory>
                {
                    FdrCategory.FullySpecific,
                    FdrCategory.SemiSpecific,
                    FdrCategory.NonSpecific
                }
            };
            Protease protease = new Protease("single N", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.Both), new Tuple<string, FragmentationTerminus>("G", FragmentationTerminus.Both) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.SingleN, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams dp = new DigestionParams(protease: protease.Name, minPeptideLength: 1, fragmentationTerminus: FragmentationTerminus.N, searchModeType: CleavageSpecificity.None);
            CommonParameters CommonParameters = new CommonParameters(
                dissociationType: DissociationType.HCD,
                precursorMassTolerance: new PpmTolerance(5),
                digestionParams: dp,
                scoreCutoff: 1,
                addCompIons: true);

            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();

            var proteinList = new List<Protein> { new Protein("GGGGGMNNNKQQQGGGGG", "TestProtein") };

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, CommonParameters, 100000, false, new List<FileInfo>(), new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndex;

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[][] allPsmsArrays = new PeptideSpectralMatch[3][];
            allPsmsArrays[0] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[1] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[2] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            PeptideSpectralMatch[] allPsmsArray = allPsmsArrays[2];
            new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, fragmentIndexDict, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray.Length);

            Assert.IsTrue(allPsmsArray[0].Score > 4);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);
            CommonParameters = new CommonParameters(
                 digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 1),
                 precursorMassTolerance: new PpmTolerance(5),
                 scoreCutoff: 1);

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            allPsmsArray[0].ResolveAllAmbiguities();
            Assert.AreEqual("QQQGGGG", allPsmsArray[0].BaseSequence);
        }

        [Test]
        public static void TestNonSpecificEnzymeSearchEngineSingleCModifications()
        {
            SearchParameters SearchParameters = new SearchParameters
            {
                SearchTarget = true,
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                LocalFdrCategories = new List<FdrCategory>
                {
                    FdrCategory.NonSpecific
                }
            };
            DigestionParams dp = new DigestionParams("singleC", minPeptideLength: 1, fragmentationTerminus: FragmentationTerminus.C, searchModeType: CleavageSpecificity.None);

            CommonParameters CommonParameters = new CommonParameters(
                dissociationType: DissociationType.HCD,
                digestionParams: dp,
                scoreCutoff: 5,
                precursorMassTolerance: new PpmTolerance(5),
                addCompIons: true);

            PeptideWithSetModifications guiltyPwsm = new PeptideWithSetModifications("DQPKLLGIETPLPKKE", null);
            var fragments = guiltyPwsm.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both);


            var myMsDataFile = new TestDataFile(guiltyPwsm.MonoisotopicMass, fragments.Select(x => x.NeutralMass.ToMz(1)).ToArray());
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Carbamidomethyl of C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            fixedModifications.Add(mod2);

            var proteinList = new List<Protein> { new Protein("GGGGGCDQPKLLGIETPLPKKEGGGGG", null) };

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.None, CommonParameters, SearchParameters.MaxFragmentSize,true, new List<FileInfo>(), new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndex;
            var precursorIndexDict = indexResults.PrecursorIndex;

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[][] allPsmsArrays = new PeptideSpectralMatch[3][];
            allPsmsArrays[0] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[1] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[2] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            PeptideSpectralMatch[] allPsmsArray = allPsmsArrays[2];
            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>());
            var searchResults = engine.Run();

            allPsmsArray[0].ResolveAllAmbiguities();
            //Check that there is no modification hanging out on the n-terminus
            Assert.AreEqual(allPsmsArray[0].FullSequence,guiltyPwsm.FullSequence);


            proteinList = new List<Protein> { new Protein("CDQPKLLGIETPLPKKEGGGGG", null) };
            guiltyPwsm = new PeptideWithSetModifications("C[Common Fixed:Carbamidomethyl on C]DQPKLLGIETPLPKKE", new Dictionary<string, Modification> { { "Carbamidomethyl on C", mod2 } });
            fragments = guiltyPwsm.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both);
            myMsDataFile = new TestDataFile(guiltyPwsm.MonoisotopicMass, fragments.Select(x => x.NeutralMass.ToMz(1)).ToArray());
            indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.None, CommonParameters, SearchParameters.MaxFragmentSize, true, new List<FileInfo>(), new List<string>());
            indexResults = (IndexingResults)indexEngine.Run();
            precursorIndexDict = indexResults.PrecursorIndex;
            peptideIndex = indexResults.PeptideIndex;
            fragmentIndexDict = indexResults.FragmentIndex;
            listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            allPsmsArrays[2] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArray = allPsmsArrays[2];
            engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>());
            searchResults = engine.Run();
            allPsmsArray[0].ResolveAllAmbiguities();
            //Check that there is a modification hanging out on the protein n-terminus
            Assert.AreEqual(allPsmsArray[0].FullSequence, guiltyPwsm.FullSequence);

            proteinList = new List<Protein> { new Protein("GGGGGCDQPKLLGIETPLPKKEGG", null) };
            indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.None, CommonParameters, SearchParameters.MaxFragmentSize, true, new List<FileInfo>(), new List<string>());
            indexResults = (IndexingResults)indexEngine.Run();
            peptideIndex = indexResults.PeptideIndex;
            fragmentIndexDict = indexResults.FragmentIndex;
            precursorIndexDict = indexResults.PrecursorIndex;
            allPsmsArrays[2] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArray = allPsmsArrays[2];
            engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>());
            searchResults = engine.Run();
            allPsmsArray[0].ResolveAllAmbiguities();
            //Check that there is a modification hanging out on the peptide n-terminus
            Assert.AreEqual(allPsmsArray[0].FullSequence, guiltyPwsm.FullSequence);
        }

        [Test]
        public static void TestNonSpecificEnzymeSearchEngineSingleNModifications()
        {
            SearchParameters SearchParameters = new SearchParameters
            {
                SearchTarget = true,
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                LocalFdrCategories = new List<FdrCategory>
                {
                    FdrCategory.NonSpecific
                }
            };
            DigestionParams dp = new DigestionParams("singleN", minPeptideLength: 1, fragmentationTerminus: FragmentationTerminus.N, searchModeType: CleavageSpecificity.None);

            CommonParameters CommonParameters = new CommonParameters(
                dissociationType: DissociationType.HCD,
                digestionParams: dp,
                scoreCutoff: 5,
                precursorMassTolerance: new PpmTolerance(5),
                addCompIons: true);

            PeptideWithSetModifications guiltyPwsm = new PeptideWithSetModifications("DQPKLLGIETPLPKKE", null);
            var fragments = guiltyPwsm.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both);


            var myMsDataFile = new TestDataFile(guiltyPwsm.MonoisotopicMass, fragments.Select(x => x.NeutralMass.ToMz(1)).ToArray());
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            ModificationMotif.TryGetMotif("C", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "Carbamidomethyl of C", _modificationType: "Common Fixed", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 57.02146372068994);
            fixedModifications.Add(mod2);

            var proteinList = new List<Protein> { new Protein("GGDQPKLLGIETPLPKKECGGGGG", null) };

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.None, CommonParameters, SearchParameters.MaxFragmentSize, true, new List<FileInfo>(), new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndex;
            var precursorIndexDict = indexResults.PrecursorIndex;

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[][] allPsmsArrays = new PeptideSpectralMatch[3][];
            allPsmsArrays[0] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[1] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[2] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            PeptideSpectralMatch[] allPsmsArray = allPsmsArrays[2];
            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>());
            var searchResults = engine.Run();

            allPsmsArray[0].ResolveAllAmbiguities();
            //Check that there is no modification hanging out on the n-terminus
            Assert.AreEqual(allPsmsArray[0].FullSequence, guiltyPwsm.FullSequence);


            proteinList = new List<Protein> { new Protein("GGGGGDQPKLLGIETPLPKKEC", null) };
            guiltyPwsm = new PeptideWithSetModifications("GGDQPKLLGIETPLPKKEC[Common Fixed:Carbamidomethyl on C]", new Dictionary<string, Modification> { { "Carbamidomethyl on C", mod2 } });
            fragments = guiltyPwsm.Fragment(CommonParameters.DissociationType, FragmentationTerminus.Both);
            myMsDataFile = new TestDataFile(guiltyPwsm.MonoisotopicMass, fragments.Select(x => x.NeutralMass.ToMz(1)).ToArray());
            indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.None, CommonParameters, SearchParameters.MaxFragmentSize, true, new List<FileInfo>(), new List<string>());
            indexResults = (IndexingResults)indexEngine.Run();
            precursorIndexDict = indexResults.PrecursorIndex;
            peptideIndex = indexResults.PeptideIndex;
            fragmentIndexDict = indexResults.FragmentIndex;
            listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            allPsmsArrays[2] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArray = allPsmsArrays[2];
            engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>());
            searchResults = engine.Run();
            allPsmsArray[0].ResolveAllAmbiguities();
            //Check that there is a modification hanging out on the protein n-terminus
            Assert.AreEqual(allPsmsArray[0].FullSequence, guiltyPwsm.FullSequence);

            proteinList = new List<Protein> { new Protein("GGDQPKLLGIETPLPKKECGGGGG", null) };
            indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.None, CommonParameters, SearchParameters.MaxFragmentSize, true, new List<FileInfo>(), new List<string>());
            indexResults = (IndexingResults)indexEngine.Run();
            peptideIndex = indexResults.PeptideIndex;
            fragmentIndexDict = indexResults.FragmentIndex;
            precursorIndexDict = indexResults.PrecursorIndex;
            allPsmsArrays[2] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArray = allPsmsArrays[2];
            engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>());
            searchResults = engine.Run();
            allPsmsArray[0].ResolveAllAmbiguities();
            //Check that there is a modification hanging out on the peptide n-terminus
            Assert.AreEqual(allPsmsArray[0].FullSequence, guiltyPwsm.FullSequence);
        }


        [Test]
        public static void TestNonSpecificEnzymeSearchEngineSingleC()
        {
            SearchParameters SearchParameters = new SearchParameters
            {
                SearchTarget = true,
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                LocalFdrCategories = new List<FdrCategory>
                {
                    FdrCategory.FullySpecific,
                    FdrCategory.SemiSpecific,
                    FdrCategory.NonSpecific
                }
            };
            Protease protease = new Protease("single C", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.Both), new Tuple<string, FragmentationTerminus>("G", FragmentationTerminus.Both) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.SingleC, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams dp = new DigestionParams(protease: protease.Name, minPeptideLength: 1, fragmentationTerminus: FragmentationTerminus.C, searchModeType: CleavageSpecificity.None);

            CommonParameters CommonParameters = new CommonParameters(
                dissociationType: DissociationType.HCD,
                digestionParams: dp,
                scoreCutoff: 4,
                precursorMassTolerance: new PpmTolerance(5),
                addCompIons: true);

            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();

            var proteinList = new List<Protein> { new Protein("GGGGGMNNNKQQQGGGGG", null) };

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, CommonParameters, SearchParameters.MaxFragmentSize, false, new List<FileInfo>(), new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndex;

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[][] allPsmsArrays = new PeptideSpectralMatch[3][];
            allPsmsArrays[0] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[1] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[2] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            PeptideSpectralMatch[] allPsmsArray = allPsmsArrays[2];
            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, fragmentIndexDict, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>());
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            //Single ms2 scan
            Assert.AreEqual(1, allPsmsArray.Length);
            Assert.That(allPsmsArray[0] != null);

            Assert.IsTrue(allPsmsArray[0].Score > 7);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);
            allPsmsArray[0].ResolveAllAmbiguities();
            Assert.AreEqual("QQQGGGG", allPsmsArray[0].BaseSequence);
        }

        [Test]
        public static void TestNonSpecificEnzymeVariableModificationHandlingNTerm()
        {
            var protein = new Protein("MGGGGGMNNNKQQQMGGGGMGM", "TestProtein");
            var protease = new Protease("singleN2", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.Both), new Tuple<string, FragmentationTerminus>("G", FragmentationTerminus.Both), new Tuple<string, FragmentationTerminus>("M", FragmentationTerminus.Both), new Tuple<string, FragmentationTerminus>("N", FragmentationTerminus.Both), new Tuple<string, FragmentationTerminus>("Q", FragmentationTerminus.Both) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.SingleN, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            var variableModifications = new List<Modification> { new Modification(_originalId: "16", _target: motifM, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915) };
            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 5, maxModsForPeptides: 3);
            var ListOfModifiedPeptides = protein.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.AreEqual(ListOfModifiedPeptides.Count, 192);

            var protein2 = new Protein(new string("MGGGGGMNNNKQQQMGGGGMGM".ToCharArray().Reverse().ToArray()), "TestProtein");
            var ListOfModifiedPeptides2 = protein2.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.AreEqual(ListOfModifiedPeptides2.Count, 132);
        }

        [Test]
        public static void TestNonSpecificEnzymeVariableModificationHandlingCTerm()
        {
            var protein = new Protein("MGGGGGMNNNKQQQMGGGGMGM", "TestProtein");
            var protease = new Protease("singleC2", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.Both), new Tuple<string, FragmentationTerminus>("G", FragmentationTerminus.Both), new Tuple<string, FragmentationTerminus>("M", FragmentationTerminus.Both), new Tuple<string, FragmentationTerminus>("N", FragmentationTerminus.Both), new Tuple<string, FragmentationTerminus>("Q", FragmentationTerminus.Both) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.SingleC, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            var variableModifications = new List<Modification> { new Modification(_originalId: "16", _target: motifM, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915) };
            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 5, maxModsForPeptides: 3);
            var ListOfModifiedPeptides = protein.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.AreEqual(ListOfModifiedPeptides.Count, 132);

            var protein2 = new Protein(new string("MGGGGGMNNNKQQQMGGGGMGM".ToCharArray().Reverse().ToArray()), "TestProtein");
            var ListOfModifiedPeptides2 = protein2.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
            Assert.AreEqual(ListOfModifiedPeptides2.Count, 192);
        }

        [Test]
        public static void TestSemiSpecificEnzymeEngineSingleN()
        {
            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();
            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
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

            var proteinList = new List<Protein> { new Protein("GGGGGMNNNKQQQGGGGGGKKRKG", "TestProtein") };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new SinglePpmAroundZeroSearchMode(5);
            var protease = new Protease("SingleN", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("K", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.None, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters(
                dissociationType: DissociationType.HCD,
                productMassTolerance: productMassTolerance,
                digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 5, maxModsForPeptides: 2, fragmentationTerminus: FragmentationTerminus.N, searchModeType: CleavageSpecificity.Semi),
                scoreCutoff: 2,
                addCompIons: true);

            HashSet<DigestionParams> digestParams = new HashSet<DigestionParams> { CommonParameters.DigestionParams };
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, CommonParameters, 100000, true, new List<FileInfo>(), new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndex;
            var precursorIndexDict = indexResults.PrecursorIndex;

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[][] allPsmsArrays = new PeptideSpectralMatch[2][];
            allPsmsArrays[0] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[1] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            PeptideSpectralMatch[] allPsmsArray = allPsmsArrays[1];
            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, 1, CommonParameters, searchModes, 0, new List<string>());
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray.Length);

            Assert.IsTrue(allPsmsArray[0].Score > 4);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);
            allPsmsArray[0].ResolveAllAmbiguities();
            Assert.AreEqual("QQQGGGG", allPsmsArray[0].BaseSequence);
        }

        [Test]
        public static void TestSemiSpecificEnzymeEngineSingleC()
        {
            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();
            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
            foreach (var mod in fixedModifications)
            {
                modsDictionary.Add(mod, 0);
            }
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

            var proteinList = new List<Protein> { new Protein("GGGGGMKNNNQQQGGGGKGG", null, null, null, null, new List<ProteolysisProduct> { new ProteolysisProduct(null, null, "test") }) };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new SinglePpmAroundZeroSearchMode(5);
            Protease protease = new Protease("SingleC", new List<Tuple<string, FragmentationTerminus>>(), new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.SingleC, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters(
                scoreCutoff: 1,
                productMassTolerance: productMassTolerance,
                digestionParams: new DigestionParams(protease: protease.Name, maxMissedCleavages: 5, minPeptideLength: 5, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C),
                dissociationType: DissociationType.HCD,
                addCompIons: true);

            HashSet<DigestionParams> digestParams = new HashSet<DigestionParams> { CommonParameters.DigestionParams };
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.Reverse, CommonParameters, 30000, false, new List<FileInfo>(), new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndex;

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[][] allPsmsArrays = new PeptideSpectralMatch[3][];
            allPsmsArrays[0] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[1] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            allPsmsArrays[2] = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            PeptideSpectralMatch[] allPsmsArray = allPsmsArrays[2];
            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArrays, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, fragmentIndexDict, 1, CommonParameters, searchModes, 0, new List<string>());
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray.Length);

            Assert.IsTrue(allPsmsArray[0].Score > 4);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);
            allPsmsArray[0].ResolveAllAmbiguities();
            Assert.AreEqual("QQQGGGG", allPsmsArray[0].BaseSequence);
        }

        [Test]
        public static void TestClassicSemiProtease()
        {
            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();
            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
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

            var proteinList = new List<Protein> { new Protein("MGGGGGMKNNNQQQGGGGKGKKNKKGN", "hello") };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new SinglePpmAroundZeroSearchMode(5);
            var protease = new Protease("semi-trypsin1", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("G", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Semi, null, null, null);
            var protease2 = new Protease("semi-trypsin2", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("N", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Semi, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);
            CommonParameters CommonParameters = new CommonParameters(
                productMassTolerance: productMassTolerance,
                digestionParams: new DigestionParams(protease: protease.Name, maxMissedCleavages: 5)
                );

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters, new List<string>()).Run();

            //////////////////////////////

            CommonParameters CommonParameters2 = new CommonParameters
            (
                productMassTolerance: productMassTolerance,
                digestionParams: new DigestionParams(protease: protease2.Name, maxMissedCleavages: 5)
            );

            HashSet<DigestionParams> digestParams2 = new HashSet<DigestionParams> { CommonParameters2.DigestionParams };

            bool DoPrecursorDeconvolution2 = true;
            bool UseProvidedPrecursorInfo2 = true;
            double DeconvolutionIntensityRatio2 = 4;
            int DeconvolutionMaxAssumedChargeState2 = 10;
            Tolerance DeconvolutionMassTolerance2 = new PpmTolerance(5);

            var listOfSortedms2Scans2 = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray2 = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray2, listOfSortedms2Scans2, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters2, new List<string>()).Run();
            // Single search mode
            Assert.AreEqual(1, allPsmsArray2.Length);
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray2.Length);
            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

            Assert.IsTrue(allPsmsArray2[0].Score > 4);
            Assert.IsTrue(allPsmsArray[0].Score > 4);
            Assert.AreEqual(2, allPsmsArray2[0].ScanNumber);
            Assert.AreEqual(allPsmsArray[0].ScanNumber, allPsmsArray2[0].ScanNumber);

            Assert.AreEqual("QQQGGGG", allPsmsArray2[0].BaseSequence);
            Assert.AreEqual(allPsmsArray[0].BaseSequence, allPsmsArray2[0].BaseSequence);
        }

        [Test]
        public static void TestClassicSemiProteolysis()
        {
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var localizeableModifications = new List<Modification>();
            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
            foreach (var mod in fixedModifications)
            {
                modsDictionary.Add(mod, 0);
            }
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
            List<ProteolysisProduct> protprod = new List<ProteolysisProduct> { new ProteolysisProduct(9, 21, "chain") };
            var proteinList = new List<Protein> { new Protein("MGGGGGMKNNNQQQGGGGKLKGKKNKKGN", "hello", null, null, null, protprod) };

            var protease = new Protease("semi-Trypsin", new List<Tuple<string, FragmentationTerminus>> { new Tuple<string, FragmentationTerminus>("G", FragmentationTerminus.C) }, new List<Tuple<string, FragmentationTerminus>>(), CleavageSpecificity.Semi, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams digestParams = new DigestionParams(protease: protease.Name, minPeptideLength: 2);

            //expect NNNQQQ, NNNQQ, NNNQ, NNN, NN and LK, KLK
            Dictionary<string, bool> found = new Dictionary<string, bool>
            {
                {"NNNQQQ", false},
                {"NNNQQ", false} ,
                {"NNNQ", false},
                {"NNN", false},
                {"NN", false},
                {"LK", false},
                {"KLK", false}
            };
            IEnumerable<PeptideWithSetModifications> PWSMs = proteinList[0].Digest(digestParams, new List<Modification>(), modsDictionary.Keys.ToList());
            foreach (PeptideWithSetModifications PWSM in PWSMs)
            {
                if (found.TryGetValue(PWSM.BaseSequence, out bool b))
                {
                    found[PWSM.BaseSequence] = true;
                }
            }
            foreach (KeyValuePair<string, bool> kvp in found)
            {
                Assert.IsTrue(kvp.Value);
            }
        }

        [Test]
        public static void TestClassicSearchOneNterminalModifiedPeptideOneScan()
        {
            //Create the modification and select the motif
            ModificationMotif.TryGetMotif("A", out ModificationMotif motif);
            Modification mod = new Modification(_originalId: "acetylation", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "Anywhere.");

            //Create the protein and add the mod
            List<Modification> modlist = new List<Modification> { mod };
            int modPositionInProtein = 1;
            var proteinListForSearch = new Protein("AAAGSAAVSGAGTPVAGPTGR", null, oneBasedModifications: new Dictionary<int, List<Modification>> { { modPositionInProtein, modlist } });
            List<Protein> proteins = new List<Protein> { proteinListForSearch };

            //Prepare the scan
            double[] mz = { 101.814659, 101.851852, 101.856758, 101.920807, 101.957397, 101.993172, 102.029053, 102.035851, 102.065628, 102.102516, 102.123955, 102.129654, 102.136391, 102.144806, 102.173798, 102.182594, 102.19574, 102.201569, 102.208534, 102.227333, 102.231323, 102.237877, 102.266922, 110.071281, 124.284775, 125.057686, 129.101776, 136.075348, 141.059647, 141.067139, 147.112045, 152.296143, 157.095963, 159.321686, 159.391266, 159.461609, 159.505859, 159.530792, 159.567841, 159.673218, 159.688232, 159.708939, 159.717117, 159.744415, 159.883804, 169.132446, 175.118637, 175.639145, 185.091217, 185.629135, 185.649612, 198.086334, 202.961533, 225.122726, 230.112839, 232.139282, 247.539169, 256.128204, 257.129761, 269.123596, 283.467743, 287.133331, 304.127319, 306.832855, 307.01059, 307.281647, 313.149384, 314.151093, 329.183044, 358.160583, 360.700592, 365.218567, 370.208008, 372.185425, 373.445496, 382.170685, 383.1745, 400.180878, 401.178406, 415.016876, 425.206726, 430.238647, 431.253052, 438.198273, 439.226074, 442.213654, 453.207336, 454.210938, 467.186005, 470.226654, 471.218048, 472.219635, 487.260254, 488.259033, 498.980499, 524.244202, 525.244385, 528.277832, 541.27124, 542.254883, 543.258118, 558.296448, 559.29895, 602.023743, 637.564087, 638.122498, 638.68103, 639.24469, 641.322876, 642.328674, 651.320923, 665.653137, 665.812622, 678.225586, 678.837524, 684.52533, 699.105225, 702.819702, 710.320374, 718.544617, 731.125793, 731.804626, 732.479858, 732.635559, 745.458618, 754.418152, 755.42041, 770.024719, 835.083557, 855.465454, 883.47406, 884.485962, 885.430664, 894.471252, 912.484192, 913.487305, 983.525391, 984.515137, 1022.541321, 1040.545166, 1041.545166, 1109.568115, 1110.565918, 1122.717896, 1127.574707, 1128.582764, 1143.576172, 1192.367432, 1209.613403, 1226.649048, 1227.651123, 1255.27356, 1270.274414, 1273.796753, 1297.668091, 1298.676636, 1520.323364, 1541.422729, 1546.209839, 1639.540283, 1669.039673, 1670.780518, 1671.45459, 1671.927368, 1672.41272, 1673.799194, 1674.881836, 1687.070557 };
            double[] intensities = { 2084.605469, 1503.049316, 1508.572144, 1997.866089, 3907.81958, 3168.024902, 2187.750488, 2006.690186, 1900.748047, 2346.96875, 5075.943359, 1605.881958, 3754.241699, 2652.021484, 2389.109375, 1921.068115, 7624.829102, 2316.127441, 4926.003906, 4061.111328, 4849.558105, 1925.951172, 2153.752441, 1949.42981, 2058.943604, 1687.470215, 1928.76416, 4824.449219, 1376.267334, 1439.787109, 2010.894897, 1483.543213, 1929.985596, 3110.926025, 3372.582275, 2304.003174, 2679.970215, 2560.099854, 1692.449585, 4134.439453, 3104.67627, 3688.150146, 3676.4646, 2814.727051, 1772.865723, 2332.302002, 4913.470703, 1477.579102, 19578.2793, 2073.369141, 2618.062988, 4346.097168, 1499.773193, 1665.242065, 2426.635742, 7056.563477, 1847.738281, 34847.85547, 2724.011963, 6279.092773, 1623.552734, 6072.998535, 3376.6521, 1946.103149, 2381.508301, 1723.251953, 43135.01953, 2463.636719, 2124.577393, 3385.583008, 1918.600098, 3191.296387, 1865.107422, 3907.073975, 2091.009766, 47135.79688, 5113.383301, 19901.54883, 1716.9729, 1691.421753, 2594.325928, 6075.093262, 1948.982544, 1560.634033, 2563.348145, 1857.768311, 96429.07813, 12289.67578, 2152.950928, 1781.416992, 49654.84766, 7577.02832, 19761.38867, 3075.865967, 2130.962158, 6758.932617, 2022.507935, 4839.808105, 2863.790039, 67360.53906, 11148.44824, 23283.18945, 4307.828125, 1818.466187, 4734.982422, 2687.210449, 4251.597168, 3365.632813, 14722.28027, 3734.265869, 1774.239624, 1827.377075, 2492.476074, 9175.493164, 2313.220215, 2536.213867, 2223.19165, 2145.397461, 2189.492676, 3675.263184, 2770.946289, 4485.819336, 2932.245605, 1713.99939, 1857.43042, 45807.97266, 10789.84961, 1848.404297, 1901.754517, 10514.50977, 14093.94141, 5819.866699, 6942.87207, 3736.169189, 34573.74219, 11260.63867, 13557.64746, 3463.163818, 7701.758301, 39518.83203, 15681.95898, 13644.61426, 7358.266602, 5733.412598, 46497.80469, 20523.03516, 1688.138916, 1865.104248, 2219.646729, 9009.290039, 4519.379395, 2556.661621, 1627.402466, 2145.759766, 4486.545898, 3421.69873, 1836.025391, 1980.848999, 2378.1521, 5128.462402, 2796.475586, 2296.582764, 8305.148438, 9054.726563, 6919.573242, 3180.789063, 2163.666504, 1787.661133 };

            MsDataScan[] msDataScanArray = new MsDataScan[1];
            msDataScanArray[0] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + 1);

            Ms2ScanWithSpecificMass[] scans = new Ms2ScanWithSpecificMass[1];
            double precursorMz = 884.4541;
            int precursorCharge = 2;
            string filePath = "path";
            scans[0] = new Ms2ScanWithSpecificMass(msDataScanArray[0], precursorMz, precursorCharge, filePath, new CommonParameters());

            //set digestion parameters and digest the protein
            DigestionParams testDigestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            //set the common parameters
            Tolerance productTolerance = new PpmTolerance(20);
            Tolerance precursorTolerance = new PpmTolerance(10);
            DissociationType d = DissociationType.CID;
            CommonParameters commonParameters = new CommonParameters(dissociationType: d, precursorMassTolerance: precursorTolerance, productMassTolerance: productTolerance, digestionParams: testDigestionParams);

            //Search the scan against the protein
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[1];
            MetaMorpheusEngineResults engineResults = new ClassicSearchEngine(allPsmsArray, scans, new List<Modification>(), new List<Modification>(), proteins, new SinglePpmAroundZeroSearchMode(20), new CommonParameters(dissociationType: DissociationType.CID), new List<string>()).Run();

            //Process the results
            List<MatchedFragmentIon> matchedFragmentIonList = allPsmsArray.SelectMany(p => p.MatchedFragmentIons).ToList();
            List<double> neutralFragmentMasses = matchedFragmentIonList.Select(m => m.NeutralTheoreticalProduct.NeutralMass).ToList();
            List<double> massToCharges = new List<double>();
            foreach (double mass in neutralFragmentMasses)
            {
                massToCharges.Add(mass.ToMz(1));
            }
            Assert.AreEqual(20, massToCharges.Count());
        }
    }
}