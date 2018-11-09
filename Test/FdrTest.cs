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
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
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
            var digested = p.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            PeptideWithSetModifications pep1 = digested[0];
            PeptideWithSetModifications pep2 = digested[1];
            PeptideWithSetModifications pep3 = digested[2];
            PeptideWithSetModifications pep4 = digested[3];

            TestDataFile t = new TestDataFile(new List<PeptideWithSetModifications> { pep1, pep2, pep3 });

            MsDataScan mzLibScan1 = t.GetOneBasedScan(2);
            Ms2ScanWithSpecificMass scan1 = new Ms2ScanWithSpecificMass(mzLibScan1, pep1.MonoisotopicMass.ToMz(1), 1, null, new CommonParameters());
            PeptideSpectralMatch psm1 = new PeptideSpectralMatch(pep1, 0, 3, 0, scan1, digestionParams, new List<MatchedFragmentIon>());

            MsDataScan mzLibScan2 = t.GetOneBasedScan(4);
            Ms2ScanWithSpecificMass scan2 = new Ms2ScanWithSpecificMass(mzLibScan2, pep2.MonoisotopicMass.ToMz(1), 1, null, new CommonParameters());
            PeptideSpectralMatch psm2 = new PeptideSpectralMatch(pep2, 1, 2, 1, scan2, digestionParams, new List<MatchedFragmentIon>());

            MsDataScan mzLibScan3 = t.GetOneBasedScan(6);
            Ms2ScanWithSpecificMass scan3 = new Ms2ScanWithSpecificMass(mzLibScan3, pep3.MonoisotopicMass.ToMz(1), 1, null, new CommonParameters());
            PeptideSpectralMatch psm3 = new PeptideSpectralMatch(pep3, 0, 1, 2, scan3, digestionParams, new List<MatchedFragmentIon>());

            psm3.AddOrReplace(pep4, 1, 1, true, new List<MatchedFragmentIon>());

            var newPsms = new List<PeptideSpectralMatch> { psm1, psm2, psm3 };
            foreach (PeptideSpectralMatch psm in newPsms)
            {
                psm.ResolveAllAmbiguities();
            }


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
            List<Modification> variableModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            List<Modification> fixedModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => CommonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif))).ToList();

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

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            //check better when using delta
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters, new List<string>()).Run();

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.None, CommonParameters, 30000, false, new List<FileInfo>(), new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();
            MassDiffAcceptor massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

            PeptideSpectralMatch[] allPsmsArrayModern = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArrayModern, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, CommonParameters, massDiffAcceptor, 0, new List<string>()).Run();

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

            listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            //check no change when using delta
            allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters, new List<string>()).Run();

            CommonParameters = new CommonParameters(useDeltaScore: true, digestionParams: new DigestionParams(minPeptideLength: 5));

            indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, 1, DecoyType.None, CommonParameters, 30000, false, new List<FileInfo>(), new List<string>());
            indexResults = (IndexingResults)indexEngine.Run();
            massDiffAcceptor = SearchTask.GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
            allPsmsArrayModern = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ModernSearchEngine(allPsmsArrayModern, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, CommonParameters, massDiffAcceptor, 0, new List<string>()).Run();

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