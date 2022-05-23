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
        [Test]
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
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"MbrAnalysisTest\MbrTest_J3.mzML"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"MbrAnalysisTest\MbrTest_K13.mzML") };
            string fastaName = @"TestData\MbrAnalysisTest\HumanFastaSlice.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMbrAnalysisOutput");

            var engine = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("ClassicSearch", classicSearch), ("PostSearchAnalysis", postSearchTask) },
                rawSlices, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
            engine.Run();

            // Not sure what's going on here
            // Still have to determine best way to write the results of MBR analysis
            string classicPath = Path.Combine(outputFolder, @"ClassicSearch\AllPSMs.psmtsv");
            var classicPsms = File.ReadAllLines(classicPath).ToList();

        }

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


            foreach (PeptideSpectralMatch psm in allPsmsArray.Where(p => p != null))
            {
                PeptideWithSetModifications pwsm = psm.BestMatchingPeptides.First().Peptide;

                List<(string FileName, CommonParameters Parameters)> fileSpecificParameters = new() { (myFile, commonParameters) };
                PeptideSpectralMatch[] peptideSpectralMatches = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

                new MiniClassicSearchEngine(pwsm, peptideSpectralMatches, listOfSortedms2Scans, variableModifications, fixedModifications, searchModes,
                    commonParameters, fileSpecificParameters, sl, new List<string>()).Run();

                Assert.AreEqual(allPsmsArray[5].BaseSequence, peptideSpectralMatches[5].BaseSequence);
                Assert.That(peptideSpectralMatches[5].SpectralAngle, Is.EqualTo(allPsmsArray[5].SpectralAngle).Within(0.01));
            }
        }
    }
}