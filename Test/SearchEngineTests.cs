using EngineLayer;
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
            CommonParameters CommonParameters = new CommonParameters
            {
                Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
                MinPeptideLength = null,
                ConserveMemory = false,
                ScoreCutoff = 1,
            };

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", null) };

            var searchModes = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5) };

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            Psm[][] allPsmsArray = new Psm[searchModes.Count()][];
            for (int aede = 0; aede < searchModes.Count; aede++)
                allPsmsArray[aede] = new Psm[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, false, CommonParameters, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // One scan
            Assert.AreEqual(1, allPsmsArray[0].Length);

            Assert.IsTrue(allPsmsArray[0][0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0][0].ScanNumber);

            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm>[] { new List<Psm> { allPsmsArray[0][0] } }, proteinList, fixedModifications, variableModifications, TerminusType.None, CommonParameters, new List<string>()).Run();

            foreach (var huh in allPsmsArray[0])
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

            Assert.AreEqual("QQQ", allPsmsArray[0][0].BaseSequence);
        }

        [Test]
        public static void TestClassicSearchEngineWithWeirdPeptide()
        {
            CommonParameters CommonParameters = new CommonParameters
            {
                Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
                MinPeptideLength = null,
                ConserveMemory = false,
                ScoreCutoff = 1,
                MaxMissedCleavages = 0,
            };

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("QXQ", null) };

            var searchModes = new List<MassDiffAcceptor> { new OpenSearchMode() };

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            Psm[][] allPsmsArray = new Psm[searchModes.Count()][];
            for (int aede = 0; aede < searchModes.Count; aede++)
                allPsmsArray[aede] = new Psm[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, false, CommonParameters, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // One Scan
            Assert.AreEqual(1, allPsmsArray[0].Length);

            Assert.IsTrue(allPsmsArray[0][0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0][0].ScanNumber);

            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm>[] { new List<Psm> { allPsmsArray[0][0] } }, proteinList, fixedModifications, variableModifications, TerminusType.None, CommonParameters, new List<string>()).Run();

            foreach (var huh in allPsmsArray[0])
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

            Assert.AreEqual("QXQ", allPsmsArray[0][0].BaseSequence);
        }

        [Test]
        public static void TestModernSearchEngine()
        {
            SearchParameters SearchParameters = new SearchParameters();
            CommonParameters CommonParameters = new CommonParameters
            {
                Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
                MinPeptideLength = null,
                ConserveMemory = false,
                ScoreCutoff = 1,
            };

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

            var searchModes = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5) };
            SearchParameters.MassDiffAcceptors = searchModes;

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B, ProductType.Y }, 1, true, CommonParameters, new List<string>());
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
            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, new List<ProductType>(), 0, CommonParameters, SearchParameters.AddCompIons, SearchParameters.MassDiffAcceptors, new List<string>()).Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray[0].Length);

            Assert.IsTrue(allPsmsArray[0][0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0][0].ScanNumber);

            var hah = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm>[] { new List<Psm> { allPsmsArray[0][0] } }, proteinList, fixedModifications, variableModifications, TerminusType.None, CommonParameters, new List<string>()).Run();

            foreach (var huh in allPsmsArray[0])
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

            Assert.AreEqual("QQQ", allPsmsArray[0][0].BaseSequence);
        }

        [Test]
        public static void TestModernSearchEngineWithWeirdPeptide()
        {
            SearchParameters SearchParameters = new SearchParameters();
            CommonParameters CommonParameters = new CommonParameters
            {
                Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
                MinPeptideLength = null,
                ConserveMemory = false,
                ScoreCutoff = 1,
            };

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var localizeableModifications = new List<ModificationWithMass>();

            var proteinList = new List<Protein> { new Protein("MNNNKQXQ", null) };

            var searchModes = new List<MassDiffAcceptor> { new OpenSearchMode() };
            SearchParameters.MassDiffAcceptors = searchModes;

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B, ProductType.Y }, 1, true, CommonParameters, new List<string>());
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
            var engine = new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, new List<ProductType>(), 0, CommonParameters, SearchParameters.AddCompIons, SearchParameters.MassDiffAcceptors, new List<string>());
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(1, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray[0].Length);

            Assert.IsTrue(allPsmsArray[0][0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0][0].ScanNumber);

            new SequencesToActualProteinPeptidesEngine(new List<Psm>[] { new List<Psm> { allPsmsArray[0][0] } }, proteinList, fixedModifications, variableModifications, TerminusType.None, CommonParameters, new List<string>()).Run();

            Assert.AreEqual(3, allPsmsArray[0][0].NumDifferentCompactPeptides);
        }

        [Test]
        public static void TestNonSpecificEnzymeEngineSingleN()
        {
            SearchParameters SearchParameters = new SearchParameters();
            SearchParameters.AddCompIons = true;
            CommonParameters CommonParameters = new CommonParameters
            {
                Protease = new Protease("singleN", new List<string> { "K, G" }, new List<string>(), TerminusType.None, CleavageSpecificity.None, null, null, null),
                ConserveMemory = false,
                ScoreCutoff = 1
            };

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

            var proteinList = new List<Protein> { new Protein("GGGGGMNNNKQQQGGGGG", "TestProtein") };

            var searchModes = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5), new OpenSearchMode() };
            SearchParameters.MassDiffAcceptors = searchModes;
            CommonParameters.MinPeptideLength = null;
            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.B }, 1, true, CommonParameters, new List<string>());
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
            CommonParameters.MinPeptideLength = 5;
            var engine = new NonSpecificEnzymeEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, new List<ProductType> { ProductType.B }, 0, CommonParameters, SearchParameters.AddCompIons, SearchParameters.MassDiffAcceptors, TerminusType.N, new List<string>());
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(2, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray[1].Length);

            Assert.IsTrue(allPsmsArray[1][0].Score > 4);
            Assert.AreEqual(2, allPsmsArray[1][0].ScanNumber);
            CommonParameters.MinPeptideLength = null;
            var hah = (SequencesToActualProteinPeptidesEngineResults)new NonSpecificEnzymeSequencesToActualPeptides(new List<Psm>[] { new List<Psm> { allPsmsArray[1][0] } }, proteinList, fixedModifications, variableModifications, TerminusType.N, CommonParameters, SearchParameters.MassDiffAcceptors, new List<string>()).Run();

            foreach (var huh in allPsmsArray[1])
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

            Assert.AreEqual("QQQGGGG", allPsmsArray[1][0].BaseSequence);
        }

        [Test]
        public static void TestNonSpecificEnzymeEngineSingleC()
        {
            SearchParameters SearchParameters = new SearchParameters();
            SearchParameters.AddCompIons = true;
            CommonParameters CommonParameters = new CommonParameters
            {
                Protease = new Protease("singleC", new List<string> { "K, G" }, new List<string>(), TerminusType.None, CleavageSpecificity.None, null, null, null),
                ConserveMemory = false,
                ScoreCutoff = 1,
            };

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
            var searchModes = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5), new OpenSearchMode() };

            SearchParameters.MassDiffAcceptors = searchModes;
            CommonParameters.MinPeptideLength = null;

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType> { ProductType.Y }, 1, true, CommonParameters, new List<string>());

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
            CommonParameters.MinPeptideLength = 5;
            var engine = new NonSpecificEnzymeEngine(allPsmsArray, listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, new List<ProductType> { ProductType.Y }, 0, CommonParameters, SearchParameters.AddCompIons, SearchParameters.MassDiffAcceptors, TerminusType.C, new List<string>());
            var searchResults = engine.Run();

            // Single search mode
            Assert.AreEqual(2, allPsmsArray.Length);

            // Single ms2 scan
            Assert.AreEqual(1, allPsmsArray[1].Length);

            Assert.IsTrue(allPsmsArray[1][0].Score > 4);
            Assert.AreEqual(2, allPsmsArray[1][0].ScanNumber);

            CommonParameters.MinPeptideLength = null;
            var hah = (SequencesToActualProteinPeptidesEngineResults)new NonSpecificEnzymeSequencesToActualPeptides(new List<Psm>[] { new List<Psm> { allPsmsArray[1][0] } }, proteinList, fixedModifications, variableModifications, TerminusType.C, CommonParameters, SearchParameters.MassDiffAcceptors, new List<string>()).Run();

            foreach (var huh in allPsmsArray[1])
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(hah.CompactPeptideToProteinPeptideMatching);

            Assert.AreEqual("QQQGGGG", allPsmsArray[1][0].BaseSequence);
        }

        [Test]
        public static void TestNonSpecificEnzymeVariableModificationHandlingNTerm()
        {
            var protein = new Protein("MGGGGGMNNNKQQQMGGGGMGM", "TestProtein");
            var protease = new Protease("singleN", new List<string> { "K, G, M, N, Q" }, new List<string>(), TerminusType.None, CleavageSpecificity.None, null, null, null);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            var variableModifications = new List<ModificationWithMass> { new ModificationWithMass("16", null, motifM, TerminusLocalization.Any, 15.994915) };
            var digestedList = protein.Digest(protease, 20, 6, null, InitiatorMethionineBehavior.Variable, variableModifications);
            foreach (var peptide in digestedList)
            {
                var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(variableModifications, 20, 3).ToList();
                var PWSM = ListOfModifiedPeptides[0];
                PeptideWithSetModifications PWSMNew = new PeptideWithSetModifications(PWSM, PWSM.OneBasedStartResidueInProtein + 3, PWSM.OneBasedEndResidueInProtein - 2);
                string PWSMSequence = PWSM.Sequence;
                string PWSMNewSequence = PWSMNew.Sequence;
                char[] PWSMNewSequenceArray = PWSMNewSequence.ToCharArray();
                for (int i = 0; i < PWSMNewSequenceArray.Count(); i++)
                {
                    if (PWSMNewSequenceArray[i] == 'M')
                    {
                        Assert.IsTrue(i != PWSMNewSequenceArray.Count() - 1);
                        Assert.IsTrue(PWSMNewSequenceArray[i + 1] == '[');
                    }
                    else if (PWSMNewSequenceArray[i] == '[')
                    {
                        Assert.IsTrue(i != 0);
                    }
                }
            }
        }

        [Test]
        public static void TestNonSpecificEnzymeVariableModificationHandlingCTerm()
        {
            var protein = new Protein("MGGGGGMNNNKQQQMGGGGMGM", "TestProtein");
            var protease = new Protease("singleC", new List<string> { "K, G, M, N, Q" }, new List<string>(), TerminusType.None, CleavageSpecificity.None, null, null, null);
            ModificationMotif.TryGetMotif("M", out ModificationMotif motifM);
            var variableModifications = new List<ModificationWithMass> { new ModificationWithMass("16", null, motifM, TerminusLocalization.Any, 15.994915, null) };
            var digestedList = protein.Digest(protease, 20, 6, null, InitiatorMethionineBehavior.Variable, variableModifications);
            foreach (var peptide in digestedList)
            {
                var ListOfModifiedPeptides = peptide.GetPeptidesWithSetModifications(variableModifications, 20, 3).ToList();
                var PWSM = ListOfModifiedPeptides[0];
                PeptideWithSetModifications PWSMNew = new PeptideWithSetModifications(PWSM, PWSM.OneBasedStartResidueInProtein + 2, PWSM.OneBasedEndResidueInProtein - 3);
                string PWSMSequence = PWSM.Sequence;
                string PWSMNewSequence = PWSMNew.Sequence;
                char[] PWSMNewSequenceArray = PWSMNewSequence.ToCharArray();
                for (int i = 0; i < PWSMNewSequenceArray.Count(); i++)
                {
                    if (PWSMNewSequenceArray[i] == 'M')
                    {
                        Assert.IsTrue(i != PWSMNewSequenceArray.Count() - 1);
                        Assert.IsTrue(PWSMNewSequenceArray[i + 1] == '[');
                    }
                    else if (PWSMNewSequenceArray[i] == '[')
                    {
                        Assert.IsTrue(i != 0);
                    }
                }
            }
        }

        #endregion Public Methods
    }
}