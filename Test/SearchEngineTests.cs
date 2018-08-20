using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
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
            Protease protease = new Protease("Customized Protease", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters
                (digestionParams: new DigestionParams(
                    protease: protease.Name,
                    minPeptideLength: 1),
                scoreCutoff: 1);

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, null, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // One scan
            Assert.AreEqual(1, allPsmsArray.Length);

            Assert.IsTrue(allPsmsArray[0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<PeptideSpectralMatch>
            { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams>
            { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>()).Run();

            foreach (var huh in allPsmsArray)
                if (huh != null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

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
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("QXQ", null) };

            var searchModes = new OpenSearchMode();

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, null, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // One Scan
            Assert.AreEqual(1, allPsmsArray.Length);

            Assert.IsTrue(allPsmsArray[0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<PeptideSpectralMatch> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>()).Run();

            foreach (var huh in allPsmsArray)
                if (huh != null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

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

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType>
            { ProductType.B, ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, SearchParameters.MaxFragmentSize, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray.Length);
            Assert.That(allPsmsArray[0] != null);

            Assert.IsTrue(allPsmsArray[0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<PeptideSpectralMatch> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>()).Run();

            foreach (var huh in allPsmsArray)
                if (huh != null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

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
                    protease: "Customized Protease",
                    minPeptideLength: 1),
                scoreCutoff: 1,
                addCompIons: true);

            var proteinList = new List<Protein> { new Protein("K", null) };

            var indexEngine = new IndexingEngine(proteinList, new List<ModificationWithMass>(), new List<ModificationWithMass>(),
                new List<ProductType> { ProductType.B, ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams },
                CommonParameters, SearchParameters.MaxFragmentSize, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(new TestDataFile(), null, true, true, 4, 10, new PpmTolerance(5)).OrderBy(b => b.PrecursorMass).ToArray();

            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

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

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B, ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, SearchParameters.MaxFragmentSize, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            //Classic
            new ClassicSearchEngine(allPsmsArray, null, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters, new List<string>()).Run();

            //Modern
            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

            //NonSpecific
            new NonSpecificEnzymeSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B }, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();
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
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();

            var proteinList = new List<Protein> { new Protein("MNNNKQXQ", null) };

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B, ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, SearchParameters.MaxFragmentSize, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            var engine = new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>());
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray.Length);
            Assert.That(allPsmsArray[0] != null);

            Assert.IsTrue(allPsmsArray[0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

            new SequencesToActualProteinPeptidesEngine(new List<PeptideSpectralMatch> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications,
                new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>()).Run();

            Assert.AreEqual(3, allPsmsArray[0].NumDifferentCompactPeptides);
        }

        [Test]
        public static void TestNonSpecificEnzymeSearchEngineSingleN()
        {
            SearchParameters SearchParameters = new SearchParameters
            {
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
            };
            Protease protease = new Protease("single N", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.None), new Tuple<string, TerminusType>("G", TerminusType.None) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.SingleN, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams dp = new DigestionParams(protease: protease.Name);
            CommonParameters CommonParameters = new CommonParameters(
                precursorMassTolerance: new PpmTolerance(5),
                digestionParams: dp,
                scoreCutoff: 1,
                addCompIons: true);

            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
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

            var proteinList = new List<Protein> { new Protein("GGGGGMNNNKQQQGGGGG", "TestProtein") };

            CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 1),
                precursorMassTolerance: new PpmTolerance(5),
                scoreCutoff: 1);
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType>
            { ProductType.B }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, 100000, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndex;

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 5),
                precursorMassTolerance: new PpmTolerance(5),
                scoreCutoff: 1,
                addCompIons: true);

            new NonSpecificEnzymeSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, fragmentIndexDict, new List<ProductType> { ProductType.B }, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>()).Run();

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

            new NonSpecificEnzymeSequencesToActualPeptides(compactPeptideToProteinPeptideMatching, new List<PeptideSpectralMatch> { allPsmsArray[0] },
                proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B }, new List<DigestionParams> { CommonParameters.DigestionParams }, massDiffAcceptor, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>()).Run();

            foreach (var huh in allPsmsArray)
                if (huh != null)
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

            Assert.AreEqual("QQQGGGG", allPsmsArray[0].BaseSequence);
        }

        [Test]
        public static void TestNonSpecificEnzymeSearchEngineSingleC()
        {
            SearchParameters SearchParameters = new SearchParameters
            {
                SearchTarget = true,
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
            };
            Protease protease = new Protease("single C", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.None), new Tuple<string, TerminusType>("G", TerminusType.None) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.SingleC, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            var dp = new DigestionParams(protease: protease.Name);

            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: dp,
                scoreCutoff: 4,
                precursorMassTolerance: new PpmTolerance(5),
                addCompIons: true);

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
            CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 1),
                precursorMassTolerance: new PpmTolerance(5),
                scoreCutoff: 4);

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType>
            { ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, SearchParameters.MaxFragmentSize, new List<string>());

            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndex;

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            CommonParameters = new CommonParameters(
               digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 5),
               precursorMassTolerance: new PpmTolerance(5),
               scoreCutoff: 4);
            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, fragmentIndexDict, new List<ProductType> { ProductType.Y }, 0, CommonParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, new List<string>());
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            //Single ms2 scan
            Assert.AreEqual(1, allPsmsArray.Length);
            Assert.That(allPsmsArray[0] != null);

            Assert.IsTrue(allPsmsArray[0].Score > 4);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

            CommonParameters = new CommonParameters(
               digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 1),
               precursorMassTolerance: new PpmTolerance(5),
               scoreCutoff: 4);

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            new NonSpecificEnzymeSequencesToActualPeptides(compactPeptideToProteinPeptideMatching, new List<PeptideSpectralMatch>
            { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, massDiffAcceptor, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>()).Run();

            foreach (var huh in allPsmsArray)
            {
                if (huh != null)
                {
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
                }
            }

            Assert.AreEqual("QQQGGGG", allPsmsArray[0].BaseSequence);
        }

        [Test]
        public static void TestNonSpecificEnzymeVariableModificationHandlingNTerm()
        {
            var protein = new Protein("MGGGGGMNNNKQQQMGGGGMGM", "TestProtein");
            var protease = new Protease("singleN2", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.None), new Tuple<string, TerminusType>("G", TerminusType.None), new Tuple<string, TerminusType>("M", TerminusType.None), new Tuple<string, TerminusType>("N", TerminusType.None), new Tuple<string, TerminusType>("Q", TerminusType.None) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.SingleN, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            var variableModifications = new List<ModificationWithMass> { new ModificationWithMass("16", null, motifM, TerminusLocalization.Any, 15.994915) };
            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 5, maxModsForPeptides: 3);
            var ListOfModifiedPeptides = protein.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).ToList();
            Assert.AreEqual(ListOfModifiedPeptides.Count, 192);

            var protein2 = new Protein(new string("MGGGGGMNNNKQQQMGGGGMGM".ToCharArray().Reverse().ToArray()), "TestProtein");
            var ListOfModifiedPeptides2 = protein2.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).ToList();
            Assert.AreEqual(ListOfModifiedPeptides2.Count, 132);
        }

        [Test]
        public static void TestNonSpecificEnzymeVariableModificationHandlingCTerm()
        {
            var protein = new Protein("MGGGGGMNNNKQQQMGGGGMGM", "TestProtein");
            var protease = new Protease("singleC2", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.None), new Tuple<string, TerminusType>("G", TerminusType.None), new Tuple<string, TerminusType>("M", TerminusType.None), new Tuple<string, TerminusType>("N", TerminusType.None), new Tuple<string, TerminusType>("Q", TerminusType.None) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.SingleC, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            var variableModifications = new List<ModificationWithMass> { new ModificationWithMass("16", null, motifM, TerminusLocalization.Any, 15.994915, null) };
            DigestionParams digestionParams = new DigestionParams(protease: protease.Name, minPeptideLength: 5, maxModsForPeptides: 3);
            var ListOfModifiedPeptides = protein.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).ToList();
            Assert.AreEqual(ListOfModifiedPeptides.Count, 132);

            var protein2 = new Protein(new string("MGGGGGMNNNKQQQMGGGGMGM".ToCharArray().Reverse().ToArray()), "TestProtein");
            var ListOfModifiedPeptides2 = protein2.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).ToList();
            Assert.AreEqual(ListOfModifiedPeptides2.Count, 192);
        }

        [Test]
        public static void TestSemiSpecificEnzymeEngineSingleN()
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

            var proteinList = new List<Protein> { new Protein("GGGGGMNNNKQQQGGGGGGKKRKG", "TestProtein") };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new SinglePpmAroundZeroSearchMode(5);
            var protease = new Protease("singleN3", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.None, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters(
                productMassTolerance: productMassTolerance,
                digestionParams: new DigestionParams(protease: protease.Name, minPeptideLength: 5, maxModsForPeptides: 2, semiProteaseDigestion: true),
                yIons: false,
                scoreCutoff: 2,
                addCompIons: true);

            HashSet<DigestionParams> digestParams = new HashSet<DigestionParams> { CommonParameters.DigestionParams };
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B }, 1, DecoyType.Reverse, digestParams, CommonParameters, 100000, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndex;

            var precursorIndexEngine = new PrecursorIndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B }, 1, DecoyType.Reverse, digestParams, CommonParameters, 100000, new List<string>());
            var precursorIndexResults = (IndexingResults)precursorIndexEngine.Run();
            var precursorIndexDict = precursorIndexResults.FragmentIndex;

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, new List<ProductType> { ProductType.B }, 1, CommonParameters, searchModes, 0, new List<string>());
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray.Length);

            Assert.IsTrue(allPsmsArray[0].Score > 4);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            new NonSpecificEnzymeSequencesToActualPeptides(compactPeptideToProteinPeptideMatching, new List<PeptideSpectralMatch> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B }, digestParams, searchModes, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>()).Run();

            foreach (var huh in allPsmsArray)
            {
                if (huh != null)
                {
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
                }
            }

            Assert.AreEqual("QQQGGGG", allPsmsArray[0].BaseSequence);
        }

        [Test]
        public static void TestSemiSpecificEnzymeEngineSingleC()
        {
            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
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
            var protease = new Protease("singleC3", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("G", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.None, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters(
                scoreCutoff: 1,
                productMassTolerance: productMassTolerance,
                digestionParams: new DigestionParams(protease: protease.Name, maxMissedCleavages: 5, minPeptideLength: 5, semiProteaseDigestion: true, terminusTypeSemiProtease: TerminusType.C),
                bIons: false,
                addCompIons: true);

            HashSet<DigestionParams> digestParams = new HashSet<DigestionParams> { CommonParameters.DigestionParams };
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.Y }, 1, DecoyType.Reverse, digestParams, CommonParameters, 30000, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();
            var peptideIndex = indexResults.PeptideIndex;
            var fragmentIndexDict = indexResults.FragmentIndex;

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, fragmentIndexDict, new List<ProductType> { ProductType.Y }, 1, CommonParameters, searchModes, 0, new List<string>());
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray.Length);

            Assert.IsTrue(allPsmsArray[0].Score > 4);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            new NonSpecificEnzymeSequencesToActualPeptides(compactPeptideToProteinPeptideMatching, new List<PeptideSpectralMatch> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.Y }, digestParams, searchModes, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>()).Run();

            foreach (var huh in allPsmsArray)
            {
                if (huh != null)
                {
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
                }
            }

            Assert.AreEqual("QQQGGGG", allPsmsArray[0].BaseSequence);
        }

        [Test]
        public static void TestClassicSemiProtease()
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

            var proteinList = new List<Protein> { new Protein("MGGGGGMKNNNQQQGGGGKGKKNKKGN", "hello") };

            var productMassTolerance = new AbsoluteTolerance(0.01);
            var searchModes = new SinglePpmAroundZeroSearchMode(5);
            var protease = new Protease("semi-trypsin1", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("G", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Semi, null, null, null);
            var protease2 = new Protease("semi-trypsin2", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("N", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Semi, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            ProteaseDictionary.Dictionary.Add(protease2.Name, protease2);
            CommonParameters CommonParameters = new CommonParameters(
                productMassTolerance: productMassTolerance,
                digestionParams: new DigestionParams(protease: protease.Name, maxMissedCleavages: 5)
                );

            HashSet<DigestionParams> digestParams = new HashSet<DigestionParams> { CommonParameters.DigestionParams };

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, null, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters, new List<string>()).Run();

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

            var listOfSortedms2Scans2 = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution2, UseProvidedPrecursorInfo2, DeconvolutionIntensityRatio2, DeconvolutionMaxAssumedChargeState2, DeconvolutionMassTolerance2).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray2 = new PeptideSpectralMatch[listOfSortedms2Scans2.Length];
            new ClassicSearchEngine(allPsmsArray2, null, listOfSortedms2Scans2, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters2, new List<string>()).Run();
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

            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<PeptideSpectralMatch> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, digestParams, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>()).Run();
            var hah2 = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<PeptideSpectralMatch> { allPsmsArray2[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, digestParams2, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>()).Run();

            foreach (var huh in allPsmsArray)
            {
                if (huh != null)
                {
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);
                }
            }
            foreach (var huh2 in allPsmsArray2)
            {
                if (huh2 != null)
                {
                    huh2.MatchToProteinLinkedPeptides(hah2.CompactPeptideToProteinPeptideMatching);
                }
            }
            Assert.AreEqual("QQQGGGG", allPsmsArray2[0].BaseSequence);
            Assert.AreEqual(allPsmsArray[0].BaseSequence, allPsmsArray2[0].BaseSequence);
        }

        [Test]
        public static void TestClassicSemiProteolysis()
        {
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
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

            var protease = new Protease("semi-Trypsin", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("G", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Semi, null, null, null);
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
            IEnumerable<PeptideWithSetModifications> PWSMs = proteinList[0].Digest(digestParams, new List<ModificationWithMass>(), modsDictionary.Keys.ToList());
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
    }
}