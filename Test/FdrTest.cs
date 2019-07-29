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
using System;
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

            psm3.AddOrReplace(pep4, 1, 1, true, new List<MatchedFragmentIon>(), 0);

            var newPsms = new List<PeptideSpectralMatch> { psm1, psm2, psm3 };
            foreach (PeptideSpectralMatch psm in newPsms)
            {
                psm.ResolveAllAmbiguities();
            }

            FdrAnalysisEngine fdr = new FdrAnalysisEngine(newPsms, searchModes.NumNotches, new CommonParameters(), nestedIds);

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

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            //check better when using delta
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, proteinList, searchModes, CommonParameters, new List<string>()).Run();

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, 1, DecoyType.None, CommonParameters, 30000, false, new List<FileInfo>(), new List<string>());
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
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, proteinList, searchModes, CommonParameters, new List<string>()).Run();

            CommonParameters = new CommonParameters(useDeltaScore: true, digestionParams: new DigestionParams(minPeptideLength: 5));

            indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, 1, DecoyType.None, CommonParameters, 30000, false, new List<FileInfo>(), new List<string>());
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

        [Test]
        public static void TestComputePEPValue()
        {
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var origDataFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams());
            var myMsDataFile = myFileManager.LoadFile(origDataFile, CommonParameters);
            var searchModes = new SinglePpmAroundZeroSearchMode(5);
            List<Protein> proteinList = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta"), true, DecoyType.Reverse, false, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                    ProteinDbLoader.UniprotOrganismRegex, out var dbErrors, -1);
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML", CommonParameters).OrderBy(b => b.PrecursorMass).ToArray();
            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, proteinList, searchModes, CommonParameters, new List<string>()).Run();
            FdrAnalysisResults fdrResultsClassicDelta = (FdrAnalysisResults)(new FdrAnalysisEngine(allPsmsArray.Where(p => p != null).ToList(), 1, CommonParameters, new List<string>()).Run());

            var nonNullPsms = allPsmsArray.Where(p => p != null).ToList();
            var nonNullPsmsOriginalCopy = allPsmsArray.Where(p => p != null).ToList();

            var maxScore = nonNullPsms.Select(n => n.Score).Max();
            var maxScorePsm = nonNullPsms.Where(n => n.Score == maxScore).First();

            Dictionary<string, int> sequenceToPsmCount = new Dictionary<string, int>();

            List<string> sequences = new List<string>();
            foreach (PeptideSpectralMatch psm in nonNullPsms)
            {
                var ss = psm.BestMatchingPeptides.Select(b => b.Peptide.FullSequence).ToList();
                sequences.Add(String.Join("|", ss));
            }

            var s = sequences.GroupBy(i => i);

            foreach (var grp in s)
            {
                sequenceToPsmCount.Add(grp.Key, grp.Count());
            }

            Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificRetTimeHI_behavior = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();
            Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificRetTemHI_behaviorModifiedPeptides = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();

            //average hydrophobicity, standard deviation hydrophobicity
            Tuple<double, double> at = new Tuple<double, double>(33.0,1.0);

            Dictionary<int, Tuple<double, double>> HI_Time_avg_dev = new Dictionary<int, Tuple<double, double>>
            {
                { 154, at }
            };

            fileSpecificRetTimeHI_behavior.Add(@"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML", HI_Time_avg_dev);


            var maxPsmData = PEP_Analysis.CreateOnePsmDataEntry(maxScorePsm, sequenceToPsmCount, fileSpecificRetTimeHI_behavior, fileSpecificRetTemHI_behaviorModifiedPeptides, null, "standard", null, !maxScorePsm.IsDecoy);
            Assert.That(maxScorePsm.PeptidesToMatchingFragments.Count, Is.EqualTo(maxPsmData.Ambiguity));
            Assert.That(maxScorePsm.DeltaScore, Is.EqualTo(maxPsmData.DeltaScore).Within(0.05));
            Assert.That((float)(maxScorePsm.Score - (int)maxScorePsm.Score), Is.EqualTo(maxPsmData.Intensity).Within(0.05));
            Assert.That(maxPsmData.HydrophobicityZScore, Is.EqualTo(5.170955).Within(0.05));
            Assert.That(maxScorePsm.BestMatchingPeptides.Select(p => p.Peptide).First().MissedCleavages, Is.EqualTo(maxPsmData.MissedCleavagesCount));
            Assert.That(maxScorePsm.BestMatchingPeptides.Select(p => p.Peptide).First().AllModsOneIsNterminus.Values.Count(), Is.EqualTo(maxPsmData.ModsCount));
            Assert.That(maxScorePsm.Notch ?? 0, Is.EqualTo(maxPsmData.Notch));
            Assert.That(maxScorePsm.PsmCount, Is.EqualTo(maxPsmData.PsmCount));
            Assert.That(maxScorePsm.ScanPrecursorCharge, Is.EqualTo(maxPsmData.ScanPrecursorCharge));

            PEP_Analysis.ComputePEPValuesForAllPSMsGeneric(nonNullPsms);

            int trueCount = 0;

            foreach (var item in allPsmsArray.Where(p => p != null))
            {
                var b = item.FdrInfo.PEP;
                if (b >= 0.5)
                {
                    trueCount++;
                }
            }

            Assert.GreaterOrEqual(32, trueCount);
        }

        [Test]
        public static void TestRemoveThisAmbiguousePeptide()
        {
            Ms2ScanWithSpecificMass scanB = new Ms2ScanWithSpecificMass(
                new MsDataScan(
                    new MzSpectrum(new double[] { }, new double[] { }, false),
                    2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null),
                100, 1, null, new CommonParameters(), null);

            PeptideSpectralMatch psm1 = new PeptideSpectralMatch(new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION", "ORGANISM"), new DigestionParams(), 1, 2, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0), 0, 10, 1, scanB, new DigestionParams(), new List<MatchedFragmentIon>(), 0);

            PeptideWithSetModifications pwsm = new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION", "ORGANISM"), new DigestionParams(), 1, 2, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            psm1.AddOrReplace(pwsm, 10, 1, true, new List<MatchedFragmentIon>(), 0);

            Assert.AreEqual(2, psm1.BestMatchingPeptides.Count());

            psm1.RemoveThisAmbiguousPeptide(1, pwsm);

            Assert.AreEqual(1, psm1.BestMatchingPeptides.Count());
        }

        [Test]
        public static void TestPEPValuesReduceAmbiguity()
        {
            //TODO
            //This is hard becuase you have to have a big enough file to create a model, enough peptides with the same exact score for the same scan, and enough differences in those ambiguous peptides that one or more get deleted after pepvalue comutation.
        }

        [Test]
        public static void TestPEPValueWorksWithStoredTrainedModel()
        {
            //TODO
            //This is hard becuase you have to have a big enough file to creat a model, enough peptides with the same exact score for the same scan, and enough differences in those ambiguous peptides that one or more get deleted after pepvalue comutation.
        }
    }
}