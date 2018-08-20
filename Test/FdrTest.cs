using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class FdrTest
    {
        [Test]
        public static void FdrTestMethod()
        {
            MassDiffAcceptor searchModes = new DotMassDiffAcceptor(null, new List<double> { 0, 1.0029 }, new PpmTolerance(5));
            List<string> nestedIds = new List<string>();

            Protein p = new Protein("MNKNNKNNNKNNNNK", null);
            DigestionParams digestionParams = new DigestionParams();
            var digested = p.Digest(digestionParams, new List<ModificationWithMass>(), new List<ModificationWithMass>()).ToList();

            PeptideWithSetModifications pep1 = digested[0];
            PeptideWithSetModifications pep2 = digested[1];
            PeptideWithSetModifications pep3 = digested[2];
            PeptideWithSetModifications pep4 = digested[3];

            TestDataFile t = new TestDataFile(new List<PeptideWithSetModifications> { pep1, pep2, pep3 });

            CompactPeptide peptide1 = new CompactPeptide(pep1, TerminusType.None);
            MsDataScan mzLibScan1 = t.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scan1 = new Ms2ScanWithSpecificMass(mzLibScan1, peptide1.MonoisotopicMassIncludingFixedMods.ToMz(1), 1, null);
            PeptideSpectralMatch psm1 = new PeptideSpectralMatch(peptide1, 0, 3, 0, scan1, digestionParams);

            CompactPeptide peptide2 = new CompactPeptide(pep2, TerminusType.None);
            MsDataScan mzLibScan2 = t.GetOneBasedScan(4);
            Ms2ScanWithSpecificMass scan2 = new Ms2ScanWithSpecificMass(mzLibScan2, peptide2.MonoisotopicMassIncludingFixedMods.ToMz(1), 1, null);
            PeptideSpectralMatch psm2 = new PeptideSpectralMatch(peptide2, 1, 2, 1, scan2, digestionParams);

            CompactPeptide peptide3 = new CompactPeptide(pep3, TerminusType.None);
            MsDataScan mzLibScan3 = t.GetOneBasedScan(6);
            Ms2ScanWithSpecificMass scan3 = new Ms2ScanWithSpecificMass(mzLibScan3, peptide3.MonoisotopicMassIncludingFixedMods.ToMz(1), 1, null);
            PeptideSpectralMatch psm3 = new PeptideSpectralMatch(peptide3, 0, 1, 2, scan3, digestionParams);

            CompactPeptide peptide4 = new CompactPeptide(pep4, TerminusType.None);
            psm3.AddOrReplace(peptide4, 1, 1, true);

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                {
                    peptide1, new HashSet<PeptideWithSetModifications>{ pep1 }
                },
                {
                    peptide2, new HashSet<PeptideWithSetModifications>{ pep2 }
                },
                {
                    peptide3, new HashSet<PeptideWithSetModifications>{ pep3 }
                },
                {
                    peptide4, new HashSet<PeptideWithSetModifications>{ pep4 }
                },
            };

            psm1.MatchToProteinLinkedPeptides(matching);
            psm2.MatchToProteinLinkedPeptides(matching);
            psm3.MatchToProteinLinkedPeptides(matching);

            var newPsms = new List<PeptideSpectralMatch> { psm1, psm2, psm3 };
            CommonParameters cp = new CommonParameters(calculateEValue: true);

            FdrAnalysisEngine fdr = new FdrAnalysisEngine(newPsms, searchModes.NumNotches, cp, nestedIds);

            fdr.Run();

            Assert.AreEqual(2, searchModes.NumNotches);
            Assert.AreEqual(0, newPsms[0].FdrInfo.CumulativeDecoyNotch);
            Assert.AreEqual(1, newPsms[0].FdrInfo.CumulativeTargetNotch);
            Assert.AreEqual(0, newPsms[1].FdrInfo.CumulativeDecoyNotch);
            Assert.AreEqual(1, newPsms[1].FdrInfo.CumulativeTargetNotch);
            Assert.AreEqual(0, newPsms[2].FdrInfo.CumulativeDecoyNotch);
            Assert.AreEqual(1, newPsms[2].FdrInfo.CumulativeTargetNotch);

            Assert.AreEqual(0, newPsms[0].FdrInfo.CumulativeDecoy);
            Assert.AreEqual(1, newPsms[0].FdrInfo.CumulativeTarget);
            Assert.AreEqual(0, newPsms[1].FdrInfo.CumulativeDecoy);
            Assert.AreEqual(2, newPsms[1].FdrInfo.CumulativeTarget);
            Assert.AreEqual(0, newPsms[2].FdrInfo.CumulativeDecoy);
            Assert.AreEqual(3, newPsms[2].FdrInfo.CumulativeTarget);
        }

        [Test]
        public static void TestDeltaValues()
        {
            CommonParameters CommonParameters = new CommonParameters(scoreCutoff: 1, useDeltaScore: true, digestionParams: new DigestionParams(minPeptideLength: 5));

            SearchParameters SearchParameters = new SearchParameters
            {
                MassDiffAcceptorType = MassDiffAcceptorType.Exact,
            };
            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();

            // Generate data for files
            Protein TargetProtein1 = new Protein("TIDEANTHE", "accession1");
            Protein TargetProtein2 = new Protein("TIDELVE", "accession2");
            Protein TargetProtein3 = new Protein("TIDENIE", "accession3");
            Protein TargetProteinLost = new Protein("PEPTIDEANTHE", "accession4");
            Protein DecoyProteinFound = new Protein("PETPLEDQGTHE", "accessiond", isDecoy: true);

            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications>
            {
                TargetProtein1.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0],
                TargetProtein2.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0],
                TargetProtein3.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0],
                DecoyProteinFound.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0]
            });

            var proteinList = new List<Protein> { TargetProtein1, TargetProtein2, TargetProtein3, TargetProteinLost, DecoyProteinFound };

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            //check better when using delta
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, null, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters, new List<string>()).Run();

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType>
            { ProductType.B, ProductType.Y }, 1, DecoyType.None, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, 30000, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
            PeptideSpectralMatch[] allPsmsArrayModern = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArrayModern, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, massDiffAcceptor, 0, new List<string>()).Run();

            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching =
                new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            if (proteinList.Any())
            {
                SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(allPsmsArray.ToList(), proteinList, fixedModifications, variableModifications, new List<ProductType>
                { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>());
                var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
                compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;
            }

            foreach (var psm in allPsmsArray)
            {
                psm.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
            }
            foreach (var psm in allPsmsArrayModern)
            {
                psm.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);
            }
            FdrAnalysisResults fdrResultsClassicDelta = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.ToList(), 1, CommonParameters, new List<string>()).Run());
            FdrAnalysisResults fdrResultsModernDelta = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArrayModern.ToList(), 1, CommonParameters, new List<string>()).Run());
            Assert.IsTrue(fdrResultsClassicDelta.PsmsWithin1PercentFdr == 3);
            Assert.IsTrue(fdrResultsModernDelta.PsmsWithin1PercentFdr == 3);

            CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5));

            //check worse when using score
            FdrAnalysisResults fdrResultsClassic = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.ToList(), 1, CommonParameters, new List<string>()).Run());
            FdrAnalysisResults fdrResultsModern = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.ToList(), 1, CommonParameters, new List<string>()).Run());
            Assert.IsTrue(fdrResultsClassic.PsmsWithin1PercentFdr == 0);
            Assert.IsTrue(fdrResultsModern.PsmsWithin1PercentFdr == 0);

            //check that when delta is bad, we used the score
            // Generate data for files
            Protein DecoyProtein1 = new Protein("TLEDAGGTHE", "accession1d", isDecoy: true);
            Protein DecoyProtein2 = new Protein("TLEDLVE", "accession2d", isDecoy: true);
            Protein DecoyProtein3 = new Protein("TLEDNIE", "accession3d", isDecoy: true);
            Protein DecoyProteinShiny = new Protein("GGGGGG", "accessionShinyd", isDecoy: true);

            myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications>
            {
                TargetProtein1.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0],
                TargetProtein2.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0],
                TargetProtein3.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0],
                DecoyProteinShiny.Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[0],
            });

            proteinList = new List<Protein>
            {
                TargetProtein1, DecoyProtein1,
                TargetProtein2, DecoyProtein2,
                TargetProtein3, DecoyProtein3,
                DecoyProteinShiny,
            };

            listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            //check no change when using delta
            allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, null, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters, new List<string>()).Run();

            CommonParameters = new CommonParameters(useDeltaScore: true, digestionParams: new DigestionParams(minPeptideLength: 5));

            indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, new List<ProductType>
            { ProductType.B, ProductType.Y }, 1, DecoyType.None, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters, 30000, new List<string>());
            indexResults = (IndexingResults)indexEngine.Run();
            massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
            allPsmsArrayModern = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArrayModern, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, new List<ProductType> { ProductType.B, ProductType.Y }, 0, CommonParameters, massDiffAcceptor, 0, new List<string>()).Run();

            var compactPeptideToProteinPeptideMatching2 = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>();
            if (proteinList.Any())
            {
                SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine2 = new SequencesToActualProteinPeptidesEngine(allPsmsArray.ToList(), proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>());
                var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine2.Run();
                compactPeptideToProteinPeptideMatching2 = res.CompactPeptideToProteinPeptideMatching;
            }

            foreach (var psm in allPsmsArray)
            {
                psm.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching2);
            }
            foreach (var psm in allPsmsArrayModern)
            {
                psm.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching2);
            }
            fdrResultsClassicDelta = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.ToList(), 1, CommonParameters, new List<string>()).Run());
            fdrResultsModernDelta = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArrayModern.ToList(), 1, CommonParameters, new List<string>()).Run());
            Assert.IsTrue(fdrResultsClassicDelta.PsmsWithin1PercentFdr == 3);
            Assert.IsTrue(fdrResultsModernDelta.PsmsWithin1PercentFdr == 3);

            CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5));

            //check no change when using score
            fdrResultsClassic = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.ToList(), 1, CommonParameters, new List<string>()).Run());
            fdrResultsModern = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArrayModern.ToList(), 1, CommonParameters, new List<string>()).Run());
            Assert.IsTrue(fdrResultsClassic.PsmsWithin1PercentFdr == 3);
            Assert.IsTrue(fdrResultsModern.PsmsWithin1PercentFdr == 3);
        }
    }
}