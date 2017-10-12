//using EngineLayer;
//using EngineLayer.ClassicSearch;
//using EngineLayer.Indexing;
//using EngineLayer.ModernSearch;
//using EngineLayer.NonSpecificEnzymeSearch;
//using MzLibUtil;
//using NUnit.Framework;
//using Proteomics;
//using System.Collections.Generic;
//using System.Linq;
//using TaskLayer;
//using UsefulProteomicsDatabases;

//namespace Test
//{
//    [TestFixture]
//    public static class SearchEngineTests
//    {
//        #region Public Methods

//        [Test]
//        public static void TestClassicSearchEngine()
//        {
//            CommonParameters CommonParameters = new CommonParameters
//            {
//                DigestionParams = new DigestionParams
//                {
//                    Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
//                    MinPeptideLength = null,
//                },
//                ConserveMemory = false,
//                ScoreCutoff = 1,
//            };

//            var myMsDataFile = new TestDataFile();
//            var variableModifications = new List<ModificationWithMass>();
//            var fixedModifications = new List<ModificationWithMass>();
//            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };

//            var searchModes = new SinglePpmAroundZeroSearchMode(5);

//            bool DoPrecursorDeconvolution = true;
//            bool UseProvidedPrecursorInfo = true;
//            double DeconvolutionIntensityRatio = 4;
//            int DeconvolutionMaxAssumedChargeState = 10;
//            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

//            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

//            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];
//            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, false, CommonParameters, new List<string>()).Run();

//            // Single search mode
//            Assert.AreEqual(1, allPsmsArray.Length);

//            // One scan
//            Assert.AreEqual(1, allPsmsArray.Length);

//            Assert.IsTrue(allPsmsArray[0].Score > 1);
//            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

//            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, new List<string>()).Run();

//            foreach (var huh in allPsmsArray)
//                if (huh != null && huh.MostProbableProteinInfo == null)
//                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

//            Assert.AreEqual("QQQ", allPsmsArray[0].BaseSequence);
//        }

//        [Test]
//        public static void TestClassicSearchEngineWithWeirdPeptide()
//        {
//            CommonParameters CommonParameters = new CommonParameters
//            {
//                DigestionParams = new DigestionParams
//                {
//                    Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
//                    MinPeptideLength = null,
//                    MaxMissedCleavages = 0,
//                },
//                ConserveMemory = false,
//                ScoreCutoff = 1,
//            };

//            var myMsDataFile = new TestDataFile();
//            var variableModifications = new List<ModificationWithMass>();
//            var fixedModifications = new List<ModificationWithMass>();
//            var proteinList = new List<Protein> { new Protein("QXQ", null) };

//            var searchModes = new OpenSearchMode();

//            bool DoPrecursorDeconvolution = true;
//            bool UseProvidedPrecursorInfo = true;
//            double DeconvolutionIntensityRatio = 4;
//            int DeconvolutionMaxAssumedChargeState = 10;
//            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

//            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

//            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];
//            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, false, CommonParameters, new List<string>()).Run();

//            // Single search mode
//            Assert.AreEqual(1, allPsmsArray.Length);

//            // One Scan
//            Assert.AreEqual(1, allPsmsArray.Length);

//            Assert.IsTrue(allPsmsArray[0].Score > 1);
//            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

//            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, new List<string>()).Run();

//            foreach (var huh in allPsmsArray)
//                if (huh != null && huh.MostProbableProteinInfo == null)
//                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

//            Assert.AreEqual("QXQ", allPsmsArray[0].BaseSequence);
//        }

//        [Test]
//        public static void TestModernSearchEngine()
//        {
//            SearchParameters SearchParameters = new SearchParameters
//            {
//                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
//            };
//            CommonParameters CommonParameters = new CommonParameters
//            {
//                PrecursorMassTolerance = new PpmTolerance(5),
//                DigestionParams = new DigestionParams
//                {
//                    Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
//                    MinPeptideLength = null,
//                },
//                ConserveMemory = false,
//                ScoreCutoff = 1,
//            };

//            var myMsDataFile = new TestDataFile();
//            var variableModifications = new List<ModificationWithMass>();
//            var fixedModifications = new List<ModificationWithMass>();
//            var localizeableModifications = new List<ModificationWithMass>();
//            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
//            foreach (var mod in fixedModifications)
//                modsDictionary.Add(mod, 0);
//            int ii = 1;
//            foreach (var mod in variableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }
//            foreach (var mod in localizeableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }

//            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };

//            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B, ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, SearchParameters.MaxFragmentSize, new List<string>());
//            var indexResults = (IndexingResults)indexEngine.Run();

//            bool DoPrecursorDeconvolution = true;
//            bool UseProvidedPrecursorInfo = true;
//            double DeconvolutionIntensityRatio = 4;
//            int DeconvolutionMaxAssumedChargeState = 10;
//            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

//            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

//            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

//            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];
//            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, SearchParameters.AddCompIons, massDiffAcceptor, new List<string>()).Run();

//            // Single search mode
//            Assert.AreEqual(1, allPsmsArray.Length);

//            // Single ms2 scan
//            Assert.AreEqual(1, allPsmsArray.Length);

//            Assert.IsTrue(allPsmsArray[0].Score > 1);
//            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

//            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, new List<string>()).Run();

//            foreach (var huh in allPsmsArray)
//                if (huh != null && huh.MostProbableProteinInfo == null)
//                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

//            Assert.AreEqual("QQQ", allPsmsArray[0].BaseSequence);
//        }

//        [Test]
//        public static void TestNonViablePSM()
//        {
//            //just check if it crashes or not.
//            SearchParameters SearchParameters = new SearchParameters();
//            CommonParameters CommonParameters = new CommonParameters
//            {
//                DigestionParams = new DigestionParams
//                {
//                    Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
//                },
//                ScoreCutoff = 255
//            };

//            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

//            var myMsDataFile = new TestDataFile(true); //empty
//            var variableModifications = new List<ModificationWithMass>();
//            var fixedModifications = new List<ModificationWithMass>();
//            var localizeableModifications = new List<ModificationWithMass>();
//            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
//            foreach (var mod in fixedModifications)
//                modsDictionary.Add(mod, 0);
//            int ii = 1;
//            foreach (var mod in variableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }
//            foreach (var mod in localizeableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }

//            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };

//            var searchModes = new SinglePpmAroundZeroSearchMode(5);

//            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B, ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, SearchParameters.MaxFragmentSize, new List<string>());
//            var indexResults = (IndexingResults)indexEngine.Run();

//            bool DoPrecursorDeconvolution = true;
//            bool UseProvidedPrecursorInfo = true;
//            double DeconvolutionIntensityRatio = 4;
//            int DeconvolutionMaxAssumedChargeState = 10;
//            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

//            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

//            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];

//            //Classic
//            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, false, CommonParameters, new List<string>()).Run();

//            //Modern
//            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, SearchParameters.AddCompIons, massDiffAcceptor, new List<string>()).Run();

//            //NonSpecific
//            new NonSpecificEnzymeSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B }, 0, CommonParameters, SearchParameters.AddCompIons, massDiffAcceptor, new List<string>()).Run();
//        }

//        [Test]
//        public static void TestModernSearchEngineWithWeirdPeptide()
//        {
//            SearchParameters SearchParameters = new SearchParameters
//            {
//                MassDiffAcceptorType = MassDiffAcceptorType.Open
//            };
//            CommonParameters CommonParameters = new CommonParameters
//            {
//                DigestionParams = new DigestionParams
//                {
//                    Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
//                    MinPeptideLength = null,
//                },
//                ConserveMemory = false,
//                ScoreCutoff = 1,
//            };

//            var myMsDataFile = new TestDataFile();
//            var variableModifications = new List<ModificationWithMass>();
//            var fixedModifications = new List<ModificationWithMass>();
//            var localizeableModifications = new List<ModificationWithMass>();

//            var proteinList = new List<Protein> { new Protein("MNNNKQXQ", null) };

//            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B, ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, SearchParameters.MaxFragmentSize, new List<string>());
//            var indexResults = (IndexingResults)indexEngine.Run();

//            bool DoPrecursorDeconvolution = true;
//            bool UseProvidedPrecursorInfo = true;
//            double DeconvolutionIntensityRatio = 4;
//            int DeconvolutionMaxAssumedChargeState = 10;
//            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

//            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
//            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

//            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];
//            var engine = new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, SearchParameters.AddCompIons, massDiffAcceptor, new List<string>());
//            var searchResults = engine.Run();

//            // Single search mode
//            Assert.AreEqual(1, allPsmsArray.Length);

//            // Single ms2 scan
//            Assert.AreEqual(1, allPsmsArray.Length);

//            Assert.IsTrue(allPsmsArray[0].Score > 1);
//            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

//            new SequencesToActualProteinPeptidesEngine(new List<Psm> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, new List<string>()).Run();

//            Assert.AreEqual(3, allPsmsArray[0].NumDifferentCompactPeptides);
//        }

//        [Test]
//        public static void TestNonSpecificEnzymeSearchEngineSingleN()
//        {
//            SearchParameters SearchParameters = new SearchParameters
//            {
//                AddCompIons = true,
//                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
//            };
//            CommonParameters CommonParameters = new CommonParameters
//            {
//                PrecursorMassTolerance = new PpmTolerance(5),
//                DigestionParams = new DigestionParams
//                {
//                    Protease = new Protease("singleN", new List<string> { "K, G" }, new List<string>(), TerminusType.None, CleavageSpecificity.SingleN, null, null, null),
//                },
//                ScoreCutoff = 1,
//                ConserveMemory = false,
//            };

//            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
//            var variableModifications = new List<ModificationWithMass>();
//            var fixedModifications = new List<ModificationWithMass>();
//            var localizeableModifications = new List<ModificationWithMass>();
//            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
//            foreach (var mod in fixedModifications)
//                modsDictionary.Add(mod, 0);
//            int ii = 1;
//            foreach (var mod in variableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }
//            foreach (var mod in localizeableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }

//            var proteinList = new List<Protein> { new Protein("GGGGGMNNNKQQQGGGGG", "TestProtein") };

//            CommonParameters.DigestionParams.MinPeptideLength = null;
//            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, 100000, new List<string>());
//            var indexResults = (IndexingResults)indexEngine.Run();
//            var peptideIndex = indexResults.PeptideIndex;
//            var fragmentIndexDict = indexResults.FragmentIndex;

//            bool DoPrecursorDeconvolution = true;
//            bool UseProvidedPrecursorInfo = true;
//            double DeconvolutionIntensityRatio = 4;
//            int DeconvolutionMaxAssumedChargeState = 10;
//            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

//            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

//            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

//            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];
//            CommonParameters.DigestionParams.MinPeptideLength = 5;
//            new NonSpecificEnzymeSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, fragmentIndexDict, new List<ProductType> { ProductType.B }, 0, CommonParameters, SearchParameters.AddCompIons, massDiffAcceptor, new List<string>()).Run();

//            // Single search mode
//            Assert.AreEqual(1, allPsmsArray.Length);

//            // Single ms2 scan
//            Assert.AreEqual(1, allPsmsArray.Length);

//            Assert.IsTrue(allPsmsArray[0].Score > 4);
//            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);
//            CommonParameters.DigestionParams.MinPeptideLength = null;

//            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
//            new NonSpecificEnzymeSequencesToActualPeptides(compactPeptideToProteinPeptideMatching, new List<Psm> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B }, new List<DigestionParams> { CommonParameters.DigestionParams }, massDiffAcceptor, CommonParameters.ReportAllAmbiguity, new List<string>()).Run();

//            foreach (var huh in allPsmsArray)
//                if (huh != null && huh.MostProbableProteinInfo == null)
//                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

//            Assert.AreEqual("QQQGGGG", allPsmsArray[0].BaseSequence);
//        }

//        [Test]
//        public static void TestNonSpecificEnzymeSearchEngineSingleC()
//        {
//            SearchParameters SearchParameters = new SearchParameters
//            {
//                AddCompIons = true,
//                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
//            };
//            CommonParameters CommonParameters = new CommonParameters
//            {
//                DigestionParams = new DigestionParams
//                {
//                    Protease = new Protease("singleC", new List<string> { "K, G" }, new List<string>(), TerminusType.None, CleavageSpecificity.SingleC, null, null, null),
//                },
//                ConserveMemory = false,
//                ScoreCutoff = 4,
//                PrecursorMassTolerance = new PpmTolerance(5),
//            };

//            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
//            var variableModifications = new List<ModificationWithMass>();
//            var fixedModifications = new List<ModificationWithMass>();
//            var localizeableModifications = new List<ModificationWithMass>();
//            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
//            foreach (var mod in fixedModifications)
//                modsDictionary.Add(mod, 0);
//            int ii = 1;
//            foreach (var mod in variableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }
//            foreach (var mod in localizeableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }

//            var proteinList = new List<Protein> { new Protein("GGGGGMNNNKQQQGGGGG", null) };
//            CommonParameters.DigestionParams.MinPeptideLength = null;

//            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.Y }, 1, DecoyType.Reverse, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, SearchParameters.MaxFragmentSize, new List<string>());

//            var indexResults = (IndexingResults)indexEngine.Run();
//            var peptideIndex = indexResults.PeptideIndex;
//            var fragmentIndexDict = indexResults.FragmentIndex;

//            bool DoPrecursorDeconvolution = true;
//            bool UseProvidedPrecursorInfo = true;
//            double DeconvolutionIntensityRatio = 4;
//            int DeconvolutionMaxAssumedChargeState = 10;
//            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

//            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
//            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

//            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];
//            CommonParameters.DigestionParams.MinPeptideLength = 5;
//            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, fragmentIndexDict, new List<ProductType> { ProductType.Y }, 0, CommonParameters, SearchParameters.AddCompIons, massDiffAcceptor, new List<string>());
//            var searchResults = engine.Run();

//            // Single search mode
//            Assert.AreEqual(1, allPsmsArray.Length);

//            // Single ms2 scan
//            Assert.AreEqual(1, allPsmsArray.Length);

//            Assert.IsTrue(allPsmsArray[0].Score > 4);
//            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

//            CommonParameters.DigestionParams.MinPeptideLength = null;

//            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
//            new NonSpecificEnzymeSequencesToActualPeptides(compactPeptideToProteinPeptideMatching, new List<Psm> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, massDiffAcceptor, CommonParameters.ReportAllAmbiguity, new List<string>()).Run();

//            foreach (var huh in allPsmsArray)
//                if (huh != null && huh.MostProbableProteinInfo == null)
//                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

//            Assert.AreEqual("QQQGGGG", allPsmsArray[0].BaseSequence);
//        }

//        [Test]
//        public static void TestNonSpecificEnzymeVariableModificationHandlingNTerm()
//        {
//            var protein = new Protein("MGGGGGMNNNKQQQMGGGGMGM", "TestProtein");
//            var protease = new Protease("singleN", new List<string> { "K, G, M, N, Q" }, new List<string>(), TerminusType.None, CleavageSpecificity.None, null, null, null);
//            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
//            var variableModifications = new List<ModificationWithMass> { new ModificationWithMass("16", null, motifM, TerminusLocalization.Any, 15.994915) };
//            DigestionParams digestionParams = new DigestionParams();
//            var digestedList = protein.Digest(digestionParams, variableModifications);
//            foreach (var peptide in digestedList)
//            {
//                var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(digestionParams, variableModifications).ToList();
//                var PWSM = ListOfModifiedPeptides[0];
//                PeptideWithSetModifications PWSMNew = new PeptideWithSetModifications(PWSM, PWSM.OneBasedStartResidueInProtein + 3, PWSM.OneBasedEndResidueInProtein - 2);
//                string PWSMSequence = PWSM.Sequence;
//                string PWSMNewSequence = PWSMNew.Sequence;
//                char[] PWSMNewSequenceArray = PWSMNewSequence.ToCharArray();
//                for (int i = 0; i < PWSMNewSequenceArray.Count(); i++)
//                {
//                    if (PWSMNewSequenceArray[i] == 'M')
//                    {
//                        Assert.IsTrue(i != PWSMNewSequenceArray.Count() - 1);
//                        Assert.IsTrue(PWSMNewSequenceArray[i + 1] == '[');
//                    }
//                    else if (PWSMNewSequenceArray[i] == '[')
//                    {
//                        Assert.IsTrue(i != 0);
//                    }
//                }
//            }
//        }

//        [Test]
//        public static void TestNonSpecificEnzymeVariableModificationHandlingCTerm()
//        {
//            var protein = new Protein("MGGGGGMNNNKQQQMGGGGMGM", "TestProtein");
//            var protease = new Protease("singleC", new List<string> { "K, G, M, N, Q" }, new List<string>(), TerminusType.None, CleavageSpecificity.None, null, null, null);
//            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
//            var variableModifications = new List<ModificationWithMass> { new ModificationWithMass("16", null, motifM, TerminusLocalization.Any, 15.994915, null) };
//            DigestionParams digestionParams = new DigestionParams();
//            var digestedList = protein.Digest(digestionParams, variableModifications);
//            foreach (var peptide in digestedList)
//            {
//                var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(digestionParams, variableModifications).ToList();
//                var PWSM = ListOfModifiedPeptides[0];
//                PeptideWithSetModifications PWSMNew = new PeptideWithSetModifications(PWSM, PWSM.OneBasedStartResidueInProtein + 2, PWSM.OneBasedEndResidueInProtein - 3);
//                string PWSMSequence = PWSM.Sequence;
//                string PWSMNewSequence = PWSMNew.Sequence;
//                char[] PWSMNewSequenceArray = PWSMNewSequence.ToCharArray();
//                for (int i = 0; i < PWSMNewSequenceArray.Count(); i++)
//                {
//                    if (PWSMNewSequenceArray[i] == 'M')
//                    {
//                        Assert.IsTrue(i != PWSMNewSequenceArray.Count() - 1);
//                        Assert.IsTrue(PWSMNewSequenceArray[i + 1] == '[');
//                    }
//                    else if (PWSMNewSequenceArray[i] == '[')
//                    {
//                        Assert.IsTrue(i != 0);
//                    }
//                }
//            }
//        }

//        [Test]
//        public static void TestSemiSpecificEnzymeEngineSingleN()
//        {
//            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
//            var variableModifications = new List<ModificationWithMass>();
//            var fixedModifications = new List<ModificationWithMass>();
//            var localizeableModifications = new List<ModificationWithMass>();
//            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
//            foreach (var mod in fixedModifications)
//                modsDictionary.Add(mod, 0);
//            int ii = 1;
//            foreach (var mod in variableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }
//            foreach (var mod in localizeableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }

//            var proteinList = new List<Protein> { new Protein("GGGGGMNNNKQQQGGGGGGKKRKG", "TestProtein") };

//            var productMassTolerance = new AbsoluteTolerance(0.01);
//            var searchModes = new SinglePpmAroundZeroSearchMode(5);
//            var protease = new Protease("singleN", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.None, null, null, null);

//            CommonParameters CommonParameters = new CommonParameters
//            {
//                ProductMassTolerance = productMassTolerance,
//                YIons = false,
//            };
//            CommonParameters.DigestionParams = new DigestionParams
//            {
//                MaxMissedCleavages = 2,
//                Protease = protease,
//                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable,
//                SemiProteaseDigestion = true,
//                TerminusTypeSemiProtease = TerminusType.N
//            };
//            HashSet<DigestionParams> digestParams = new HashSet<DigestionParams> { CommonParameters.DigestionParams };
//            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B }, 1, DecoyType.Reverse, digestParams, CommonParameters, 100000, new List<string>());
//            var indexResults = (IndexingResults)indexEngine.Run();
//            var peptideIndex = indexResults.PeptideIndex;
//            var fragmentIndexDict = indexResults.FragmentIndex;

//            var precursorIndexEngine = new PrecursorIndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B }, 1, DecoyType.Reverse, digestParams, CommonParameters, 100000, new List<string>());
//            var precursorIndexResults = (IndexingResults)precursorIndexEngine.Run();
//            var precursorIndexDict = precursorIndexResults.FragmentIndex;

//            bool DoPrecursorDeconvolution = true;
//            bool UseProvidedPrecursorInfo = true;
//            double DeconvolutionIntensityRatio = 4;
//            int DeconvolutionMaxAssumedChargeState = 10;
//            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

//            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

//            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];
//            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, precursorIndexDict, new List<ProductType> { ProductType.B }, 1, CommonParameters, true, searchModes, new List<string>());
//            var searchResults = engine.Run();

//            // Single search mode
//            Assert.AreEqual(1, allPsmsArray.Length);

//            // Single ms2 scan
//            Assert.AreEqual(1, allPsmsArray.Length);

//            Assert.IsTrue(allPsmsArray[0].Score > 4);
//            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

//            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
//            new NonSpecificEnzymeSequencesToActualPeptides(compactPeptideToProteinPeptideMatching, new List<Psm> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B }, digestParams, searchModes, CommonParameters.ReportAllAmbiguity, new List<string>()).Run();

//            foreach (var huh in allPsmsArray)
//                if (huh != null && huh.MostProbableProteinInfo == null)
//                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

//            Assert.AreEqual("QQQGGGG", allPsmsArray[0].BaseSequence);
//        }

//        [Test]
//        public static void TestSemiSpecificEnzymeEngineSingleC()
//        {
//            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
//            var variableModifications = new List<ModificationWithMass>();
//            var fixedModifications = new List<ModificationWithMass>();
//            var localizeableModifications = new List<ModificationWithMass>();
//            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
//            foreach (var mod in fixedModifications)
//                modsDictionary.Add(mod, 0);
//            int ii = 1;
//            foreach (var mod in variableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }
//            foreach (var mod in localizeableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }

//            var proteinList = new List<Protein> { new Protein("GGGGGMKNNNQQQGGGGKGG", null) };

//            var productMassTolerance = new AbsoluteTolerance(0.01);
//            var searchModes = new SinglePpmAroundZeroSearchMode(5);
//            var protease = new Protease("singleC", new List<string> { "G" }, new List<string>(), TerminusType.C, CleavageSpecificity.None, null, null, null);

//            CommonParameters CommonParameters = new CommonParameters
//            {
//                ScoreCutoff = 1,
//                ProductMassTolerance = productMassTolerance,
//                BIons = false
//            };
//            CommonParameters.DigestionParams = new DigestionParams
//            {
//                MaxMissedCleavages = 5,
//                Protease = protease,
//                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable,
//                SemiProteaseDigestion = true,
//                TerminusTypeSemiProtease = TerminusType.C
//            };
//            HashSet<DigestionParams> digestParams = new HashSet<DigestionParams> { CommonParameters.DigestionParams };
//            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.Y }, 1, DecoyType.Reverse, digestParams, CommonParameters, 30000, new List<string>());
//            var indexResults = (IndexingResults)indexEngine.Run();
//            var peptideIndex = indexResults.PeptideIndex;
//            var fragmentIndexDict = indexResults.FragmentIndex;

//            bool DoPrecursorDeconvolution = true;
//            bool UseProvidedPrecursorInfo = true;
//            double DeconvolutionIntensityRatio = 4;
//            int DeconvolutionMaxAssumedChargeState = 10;
//            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

//            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

//            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];
//            var engine = new NonSpecificEnzymeSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, fragmentIndexDict, fragmentIndexDict, new List<ProductType> { ProductType.Y }, 1, CommonParameters, true, searchModes, new List<string>());
//            var searchResults = engine.Run();

//            // Single search mode
//            Assert.AreEqual(1, allPsmsArray.Length);

//            // Single ms2 scan
//            Assert.AreEqual(1, allPsmsArray.Length);

//            Assert.IsTrue(allPsmsArray[0].Score > 4);
//            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

//            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
//            new NonSpecificEnzymeSequencesToActualPeptides(compactPeptideToProteinPeptideMatching, new List<Psm> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.Y }, digestParams, searchModes, CommonParameters.ReportAllAmbiguity, new List<string>()).Run();

//            foreach (var huh in allPsmsArray)
//                if (huh != null && huh.MostProbableProteinInfo == null)
//                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

//            Assert.AreEqual("QQQGGGG", allPsmsArray[0].BaseSequence);
//        }

//        [Test]
//        public static void TestClassicSemiProtease()
//        {
//            var myMsDataFile = new TestDataFile("Yes, I'd like one slightly larger please");
//            var variableModifications = new List<ModificationWithMass>();
//            var fixedModifications = new List<ModificationWithMass>();
//            var localizeableModifications = new List<ModificationWithMass>();
//            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
//            foreach (var mod in fixedModifications)
//                modsDictionary.Add(mod, 0);
//            int ii = 1;
//            foreach (var mod in variableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }
//            foreach (var mod in localizeableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }

//            var proteinList = new List<Protein> { new Protein("MGGGGGMKNNNQQQGGGGKGKKNKKGN", "hello") };

//            var productMassTolerance = new AbsoluteTolerance(0.01);
//            var searchModes = new SinglePpmAroundZeroSearchMode(5);
//            var protease = new Protease("semi-trypsin", new List<string> { "G" }, new List<string>(), TerminusType.C, CleavageSpecificity.Semi, null, null, null);
//            var protease2 = new Protease("semi-trypsin", new List<string> { "N" }, new List<string>(), TerminusType.C, CleavageSpecificity.Semi, null, null, null);

//            CommonParameters CommonParameters = new CommonParameters
//            {
//                ProductMassTolerance = productMassTolerance,
//            };
//            CommonParameters.DigestionParams = new DigestionParams
//            {
//                MaxMissedCleavages = 5,
//                Protease = protease,
//                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable
//            };
//            HashSet<DigestionParams> digestParams = new HashSet<DigestionParams> { CommonParameters.DigestionParams };

//            bool DoPrecursorDeconvolution = true;
//            bool UseProvidedPrecursorInfo = true;
//            double DeconvolutionIntensityRatio = 4;
//            int DeconvolutionMaxAssumedChargeState = 10;
//            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

//            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

//            Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];
//            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, false, CommonParameters, new List<string>()).Run();

//            //////////////////////////////

//            CommonParameters CommonParameters2 = new CommonParameters
//            {
//                ProductMassTolerance = productMassTolerance,
//            };
//            CommonParameters2.DigestionParams = new DigestionParams
//            {
//                MaxMissedCleavages = 5,
//                Protease = protease2,
//                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Cleave
//            };
//            HashSet<DigestionParams> digestParams2 = new HashSet<DigestionParams> { CommonParameters2.DigestionParams };

//            bool DoPrecursorDeconvolution2 = true;
//            bool UseProvidedPrecursorInfo2 = true;
//            double DeconvolutionIntensityRatio2 = 4;
//            int DeconvolutionMaxAssumedChargeState2 = 10;
//            Tolerance DeconvolutionMassTolerance2 = new PpmTolerance(5);

//            var listOfSortedms2Scans2 = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution2, UseProvidedPrecursorInfo2, DeconvolutionIntensityRatio2, DeconvolutionMaxAssumedChargeState2, DeconvolutionMassTolerance2).OrderBy(b => b.PrecursorMass).ToArray();

//            Psm[] allPsmsArray2 = new Psm[listOfSortedms2Scans2.Length];
//            new ClassicSearchEngine(allPsmsArray2, listOfSortedms2Scans2, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, false, CommonParameters2, new List<string>()).Run();
//            // Single search mode
//            Assert.AreEqual(1, allPsmsArray2.Length);
//            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

//            // Single ms2 scan
//            Assert.AreEqual(1, allPsmsArray2.Length);
//            Assert.AreEqual(allPsmsArray.Length, allPsmsArray2.Length);

//            Assert.IsTrue(allPsmsArray2[0].Score > 4);
//            Assert.IsTrue(allPsmsArray[0].Score > 4);
//            Assert.AreEqual(2, allPsmsArray2[0].ScanNumber);
//            Assert.AreEqual(allPsmsArray[0].ScanNumber, allPsmsArray2[0].ScanNumber);

//            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm> { allPsmsArray[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, digestParams, CommonParameters.ReportAllAmbiguity, new List<string>()).Run();
//            var hah2 = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm> { allPsmsArray2[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, digestParams2, CommonParameters.ReportAllAmbiguity, new List<string>()).Run();

//            foreach (var huh in allPsmsArray)
//                if (huh != null && huh.MostProbableProteinInfo == null)
//                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);
//            foreach (var huh2 in allPsmsArray2)
//                if (huh2 != null && huh2.MostProbableProteinInfo == null)
//                    huh2.MatchToProteinLinkedPeptides(hah2.CompactPeptideToProteinPeptideMatching);

//            Assert.AreEqual("QQQGGGG", allPsmsArray2[0].BaseSequence);
//            Assert.AreEqual(allPsmsArray[0].BaseSequence, allPsmsArray2[0].BaseSequence);
//        }

//        [Test]
//        public static void TestClassicSemiProteolysis()
//        {
//            var variableModifications = new List<ModificationWithMass>();
//            var fixedModifications = new List<ModificationWithMass>();
//            var localizeableModifications = new List<ModificationWithMass>();
//            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
//            foreach (var mod in fixedModifications)
//                modsDictionary.Add(mod, 0);
//            int ii = 1;
//            foreach (var mod in variableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }
//            foreach (var mod in localizeableModifications)
//            {
//                modsDictionary.Add(mod, (ushort)ii);
//                ii++;
//            }
//            List<ProteolysisProduct> protprod = new List<ProteolysisProduct> { new ProteolysisProduct(9, 21, "chain") };
//            var proteinList = new List<Protein> { new Protein("MGGGGGMKNNNQQQGGGGKLKGKKNKKGN", "hello", null, null, protprod) };

//            var protease = new Protease("semi-trypsin", new List<string> { "G" }, new List<string>(), TerminusType.C, CleavageSpecificity.Semi, null, null, null);

//            DigestionParams digestParams = new DigestionParams
//            {
//                MaxMissedCleavages = 2,
//                Protease = protease,
//                MinPeptideLength = 2
//            };

//            //expect NNNQQQ, NNNQQ, NNNQ, NNN, NN and LK, KLK
//            Dictionary<string, bool> found = new Dictionary<string, bool>
//            {
//                {"NNNQQQ", false},
//                {"NNNQQ", false} ,
//                {"NNNQ", false},
//                {"NNN", false},
//                {"NN", false},
//                {"LK", false},
//                {"KLK", false}
//            };
//            IEnumerable<PeptideWithPossibleModifications> PWSMs = proteinList[0].Digest(digestParams, modsDictionary.Keys);
//            foreach (PeptideWithPossibleModifications PWSM in PWSMs)
//            {
//                if (found.TryGetValue(PWSM.BaseSequence, out bool b))
//                    found[PWSM.BaseSequence] = true;
//            }
//            foreach (KeyValuePair<string, bool> kvp in found)
//                Assert.IsTrue(kvp.Value);
//        }

//        #endregion Public Methods
//    }
//}