using Nett;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using EngineLayer;
using FlashLFQ;
using MathNet.Numerics.Statistics; // Necessary for calculating correlation 
using System.Text.RegularExpressions;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using EngineLayer.ClassicSearch;


namespace Test
{
    [TestFixture]
    public static class
    MbrAnalysisTest
    {
        /*[Test]
        public static void MbrPostSearchAnalysisTest()
        {
            SearchTask classicSearch = new SearchTask()
            {
                SearchParameters = new SearchParameters()
                {
                    MatchBetweenRuns = true,
                    WriteSpectralLibrary = true
                }
            };

            PostSearchAnalysisTask postSearchTask = new PostSearchAnalysisTask()
            {
                Parameters = new PostSearchAnalysisParameters()
                {
                    SearchParameters = new SearchParameters()
                    {
                        MatchBetweenRuns = true,
                        WriteSpectralLibrary = true
                    }
                },
                CommonParameters = new CommonParameters()
            };
            List<int> counts = new List<int>();
            List<string> rawSlices = new List<string> { 
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"MbrTestData\f1r1_sliced_mbr.raw"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"MbrTestData\f1r2_sliced_mbr.raw") };
            string fastaName = @"TestData\MbrTestData\MbrDataPrunedDB.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMbrAnalysisOutput");

            var engine = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("ClassicSearch", classicSearch), ("PostSearchAnalysis", postSearchTask) },
                rawSlices, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
            engine.Run();

            // Not sure what's going on here
            // Still have to determine best way to write the results of MBR analysis
            string classicPath = Path.Combine(outputFolder, @"ClassicSearch\AllPSMs.psmtsv");
            var classicPsms = File.ReadAllLines(classicPath).ToList();

        }*/

        [Test]
        public static void MiniClassicSearchEngineTest()
        {
            var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch");
            string myFile = Path.Combine(testDir, @"slicedMouse.raw");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters commonParameters = new CommonParameters(maxThreadsToUsePerFile: 1, scoreCutoff: 1);
            MsDataFile myMsDataFile = myFileManager.LoadFile(myFile, commonParameters);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("VIHDNFGIVEGLMTTVHAITATQK", "P16858") };

            string targetSpectralLibrary = Path.Combine(testDir, @"P16858_target.msp");
            string decoySpectralLibrary = Path.Combine(testDir, @"P16858_decoy.msp");

            List<string> specLibs = new List<string> { targetSpectralLibrary, decoySpectralLibrary };

            SpectralLibrary sl = new(specLibs);

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            Ms2ScanWithSpecificMass[] listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpectralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, commonParameters, null, sl, new List<string>(), writeSpectralLibrary).Run();

            // Single search mode
            Assert.AreEqual(7, allPsmsArray.Length);
            Assert.IsTrue(allPsmsArray[5].Score > 38);
            Assert.AreEqual("VIHDNFGIVEGLMTTVHAITATQK", allPsmsArray[5].BaseSequence);
            Assert.IsTrue(!allPsmsArray[5].IsDecoy);

            SpectralLibrarySearchFunction.CalculateSpectralAngles(sl, allPsmsArray, listOfSortedms2Scans, commonParameters);
            Assert.That(allPsmsArray[5].SpectralAngle, Is.EqualTo(0.82).Within(0.01));


            foreach (PeptideSpectralMatch psm in allPsmsArray.Where(p=>p != null))
            {
                PeptideWithSetModifications pwsm = psm.BestMatchingPeptides.First().Peptide;

                List<(string FileName, CommonParameters Parameters)> fileSpecificParameters = new() {(myFile,commonParameters) };
                PeptideSpectralMatch[] peptideSpectralMatches = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

                new MiniClassicSearchEngine(pwsm, peptideSpectralMatches, listOfSortedms2Scans, variableModifications, fixedModifications, searchModes,
                    commonParameters, fileSpecificParameters, sl, new List<string>()).Run();

                Assert.AreEqual(allPsmsArray[5].BaseSequence, peptideSpectralMatches[5].BaseSequence);
                Assert.That(peptideSpectralMatches[5].SpectralAngle,Is.EqualTo(allPsmsArray[5].SpectralAngle).Within(0.01));
            }
        }

        /*[Test]
        public static void FalseMbrMCSETest()
        {
            SearchTask classicSearch = new SearchTask()
            {
                SearchParameters = new SearchParameters()
                {
                    MatchBetweenRuns = true
                }
            };

            List<string> rawSlices = new List<string> { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"MbrTestData\f1r1_sliced_mbr.raw") };
            string fastaName = @"TestData\MbrTestData\MbrDataPrunedDB.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMbrAnalysisOutput");

            FileSpecificParameters[] fileSpecificParameters = new FileSpecificParameters[1];
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\FileSpecificActual.toml");
            var fileSpecificToml = Toml.ReadFile(filePath, MetaMorpheusTask.tomlConfig);
            fileSpecificParameters[0] = new FileSpecificParameters(fileSpecificToml);

            //put this back when it gets fixed.
            //List<PeptideSpectralMatch> allPsms = classicSearch.RunForPTMs(outputFolder, new List<DbForTask> { new DbForTask(fastaName, false) }, rawSlices, "PtmReturn", fileSpecificParameters);

            var engine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("ClassicSearch", classicSearch) }, rawSlices, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
            engine.Run();

            string classicPath = Path.Combine(outputFolder, @"ClassicSearch\AllPSMs.psmtsv");
            var classicPsms = File.ReadAllLines(classicPath).ToList();

        }*/


        /*
        public static void RealDataMbrTest()
        {
            string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"PSMsForMbrTest.psmtsv");
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"f1r1_sliced_mbr.raw"), "a", 0, 0, 0);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"f1r2_sliced_mbr.raw"), "a", 1, 0, 0);
            List<Identification> ids = new List<Identification>();
            Dictionary<string, FlashLFQ.ProteinGroup> allProteinGroups = new Dictionary<string, FlashLFQ.ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }
                SpectraFileInfo file = null;
                if (split[0].Contains("f1r1"))
                {
                    file = f1r1;
                }
                else if (split[0].Contains("f1r2"))
                {
                    file = f1r2;
                }
                string baseSequence = split[12];
                string fullSequence = split[13];
                double monoMass = double.Parse(split[21]);
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[24].Split(new char[] { '|' });
                List<FlashLFQ.ProteinGroup> proteinGroups = new List<FlashLFQ.ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new FlashLFQ.ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }
                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);
            }
            var engine = new FlashLfqEngine(ids, matchBetweenRuns: true, requireMsmsIdInCondition: false, maxThreads: 1);
            var results = engine.Run();
            var f1r1MbrResults = results
                .PeptideModifiedSequences
                .Where(p => p.Value.GetDetectionType(f1r1) == DetectionType.MBR && p.Value.GetDetectionType(f1r2) == DetectionType.MSMS).ToList();
            Assert.That(f1r1MbrResults.Count >= 132);
            var f1r2MbrResults = results.PeptideModifiedSequences
                .Where(p => p.Value.GetDetectionType(f1r1) == DetectionType.MSMS && p.Value.GetDetectionType(f1r2) == DetectionType.MBR).ToList();
            Assert.That(f1r2MbrResults.Count >= 77);
            List<(double, double)> peptideIntensities = new List<(double, double)>();
            foreach (var peptide in f1r1MbrResults)
            {
                double mbrIntensity = Math.Log(peptide.Value.GetIntensity(f1r1));
                double msmsIntensity = Math.Log(peptide.Value.GetIntensity(f1r2));
                peptideIntensities.Add((mbrIntensity, msmsIntensity));
            }
            double corr = Correlation.Pearson(peptideIntensities.Select(p => p.Item1), peptideIntensities.Select(p => p.Item2));
            Assert.That(corr > 0.8);
            peptideIntensities.Clear();
            foreach (var peptide in f1r2MbrResults)
            {
                double mbrIntensity = Math.Log(peptide.Value.GetIntensity(f1r2));
                double msmsIntensity = Math.Log(peptide.Value.GetIntensity(f1r1));
                peptideIntensities.Add((mbrIntensity, msmsIntensity));
            }
            corr = Correlation.Pearson(peptideIntensities.Select(p => p.Item1), peptideIntensities.Select(p => p.Item2));
            Assert.That(corr > 0.7);
            // the "requireMsmsIdInCondition" field requires that at least one MS/MS identification from a protein
            // has to be observed in a condition for match-between-runs
            f1r1.Condition = "b";
            engine = new FlashLfqEngine(ids, matchBetweenRuns: true, requireMsmsIdInCondition: true, maxThreads: 1);
            results = engine.Run();
            var proteinsObservedInF1 = ids.Where(p => p.FileInfo == f1r1).SelectMany(p => p.ProteinGroups).Distinct().ToList();
            var proteinsObservedInF2 = ids.Where(p => p.FileInfo == f1r2).SelectMany(p => p.ProteinGroups).Distinct().ToList();
            var proteinsObservedInF1ButNotF2 = proteinsObservedInF1.Except(proteinsObservedInF2).ToList();
            foreach (ProteinGroup protein in proteinsObservedInF1ButNotF2)
            {
                Assert.That(results.ProteinGroups[protein.ProteinGroupName].GetIntensity(f1r2) == 0);
            }
        }
        */

    }
}