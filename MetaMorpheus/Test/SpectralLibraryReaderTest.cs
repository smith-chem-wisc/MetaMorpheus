using NUnit.Framework;
using System.IO;
using System;
using System.Linq;
using EngineLayer;
using TaskLayer;
using System.Collections.Generic;
using Omics.Fragmentation;
using Proteomics;
using MassSpectrometry;
using EngineLayer.ClassicSearch;
using EngineLayer.DatabaseLoading;
using Omics.Modifications;
using Readers.SpectralLibrary;
using MzLibUtil;
using Proteomics.ProteolyticDigestion;

namespace Test
{
    [TestFixture]
    public static class TestSpectralLibrarySearch
    {
        [Test]
        public static void SpectralLibrarySearchTest()
        {
            var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch");
            var outputDir = Path.Combine(testDir, @"SpectralLibrarySearchTest");

            string library1 = Path.Combine(testDir, @"P16858_target.msp");
            string library2 = Path.Combine(testDir, @"P16858_decoy.msp");
            string fastaDb = Path.Combine(testDir, @"P16858.fasta");
            string spectraFile = Path.Combine(testDir, @"slicedMouse.raw");

            Directory.CreateDirectory(outputDir);

            var searchTask = new SearchTask();

            searchTask.RunTask(outputDir,
                new List<DbForTask>
                {
                    new DbForTask(library1, false),
                    new DbForTask(library2, false),
                    new DbForTask(fastaDb, false)
                },
                new List<string> { spectraFile },
                "");

            var results = File.ReadAllLines(Path.Combine(outputDir, @"AllPSMs.psmtsv"));
            var split = results[0].Split('\t');
            int ind = Array.IndexOf(split, "Normalized Spectral Angle");
            Assert.That(ind >= 0);

            var spectralAngle = double.Parse(results[1].Split('\t')[ind]);
            Assert.That(Math.Round(spectralAngle, 2) == 0.82);

            Directory.Delete(outputDir, true);
        }
        [Test]
        public static void DaltonToleranceSpectralLibrarySearchTest()
        {
            var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch");
            var outputDir = Path.Combine(testDir, @"DaltonToleranceSpectralLibrarySearchTest");

            string library1 = Path.Combine(testDir, @"P16858_target.msp");
            string library2 = Path.Combine(testDir, @"P16858_decoy.msp");
            string fastaDb = Path.Combine(testDir, @"P16858.fasta");
            string spectraFile = Path.Combine(testDir, @"slicedMouse.raw");

            Directory.CreateDirectory(outputDir);
            Tolerance t1 = Tolerance.ParseToleranceString("0.5 Absolute");     // AbsoluteTolerance(0.5)
            CommonParameters cp = new CommonParameters(productMassTolerance: t1);
            var searchTask = new SearchTask();
            searchTask.CommonParameters = cp;
            searchTask.RunTask(outputDir,
                new List<DbForTask>
                {
                    new DbForTask(library1, false),
                    new DbForTask(library2, false),
                    new DbForTask(fastaDb, false)
                },
                new List<string> { spectraFile },
                "");

            var results = File.ReadAllLines(Path.Combine(outputDir, @"AllPSMs.psmtsv"));
            var split = results[0].Split('\t');
            int ind = Array.IndexOf(split, "Normalized Spectral Angle");
            Assert.That(ind >= 0);

            var spectralAngle = double.Parse(results[1].Split('\t')[ind]);
            Assert.That(Math.Round(spectralAngle, 3) == 0.812);

            Directory.Delete(outputDir, true);
        }
        /// <summary>
        /// Test ensures peptide FDR is calculated and that it doesn't output PSM FDR results
        /// </summary>
        [Test]
        public static void AnotherSpectralLibrarySearchTest()
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

            SpectralLibrary sl = new SpectralLibrary(specLibs);

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpectralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, commonParameters,null, sl, new List<string>(), writeSpectralLibrary).Run();

            // Single search mode
            Assert.That(allPsmsArray.Length, Is.EqualTo(7));
            Assert.That(allPsmsArray[5].Score > 38);
            Assert.That(allPsmsArray[5].BaseSequence, Is.EqualTo("VIHDNFGIVEGLMTTVHAITATQK"));


            SpectralLibrarySearchFunction.CalculateSpectralAngles(sl, allPsmsArray, listOfSortedms2Scans, commonParameters);
            Assert.That(allPsmsArray[5].SpectralAngle, Is.EqualTo(0.82).Within(0.01));
        }


        /// <summary>
        /// Test ensures peptide FDR is calculated and that it doesn't output PSM FDR results
        /// </summary>
        [Test]
        public static void AnotherSpectralLibrarySearchTestDecoy()
        {
            var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch");
            string myFile = Path.Combine(testDir, @"slicedMouse.raw");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters commonParameters = new CommonParameters(maxThreadsToUsePerFile: 1, scoreCutoff: 1);
            MsDataFile myMsDataFile = myFileManager.LoadFile(myFile, commonParameters);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("QTATIAHVTTMLGEVIGFNDHIVK", "P16858") };

            string targetSpectralLibrary = Path.Combine(testDir, @"P16858_target.msp");
            string decoySpectralLibrary = Path.Combine(testDir, @"P16858_decoy.msp");

            List<string> specLibs = new List<string> { targetSpectralLibrary, decoySpectralLibrary };

            SpectralLibrary sl = new SpectralLibrary(specLibs);

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpectralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, commonParameters, null, sl, new List<string>(), writeSpectralLibrary).Run();

            // Single search mode
            Assert.That(allPsmsArray.Length, Is.EqualTo(7));
            Assert.That(allPsmsArray[5].Score > 38);
            Assert.That(allPsmsArray[5].BaseSequence, Is.EqualTo("VIHDNFGIVEGLMTTVHAITATQK"));
            Assert.That(allPsmsArray[5].IsDecoy);

            SpectralLibrarySearchFunction.CalculateSpectralAngles(sl, allPsmsArray, listOfSortedms2Scans, commonParameters);
            Assert.That(allPsmsArray[5].SpectralAngle, Is.EqualTo(0.82).Within(0.01));

        }
        /// <summary>
        /// Test ensures spectral angle calculation works correctly for decoys when no decoy library spectrum is available (fallback to generated decoy spectrum)
        /// </summary>
        [Test]
        public static void AnotherSpectralLibrarySearchTestDecoyNoLibrarySpectrum()
        {
            var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch");
            string myFile = Path.Combine(testDir, @"slicedMouse.raw");
            MyFileManager myFileManager = new MyFileManager(true);
            CommonParameters commonParameters = new CommonParameters(maxThreadsToUsePerFile: 1, scoreCutoff: 1);
            MsDataFile myMsDataFile = myFileManager.LoadFile(myFile, commonParameters);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("QTATIAHVTTMLGEVIGFNDHIVK", "P16858") };

            string targetSpectralLibrary = Path.Combine(testDir, @"P16858_decoy.msp"); // Using P16858_decoy.msp file which actually contains the target spectrum for this test case

            List<string> specLibs = new List<string> { targetSpectralLibrary };

            SpectralLibrary sl = new SpectralLibrary(specLibs);

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            bool writeSpectralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, commonParameters, null, sl, new List<string>(), writeSpectralLibrary).Run();

            // Single search mode
            Assert.That(allPsmsArray.Length, Is.EqualTo(7));
            Assert.That(allPsmsArray[5].Score > 38);
            Assert.That(allPsmsArray[5].BaseSequence, Is.EqualTo("VIHDNFGIVEGLMTTVHAITATQK")); //this is the decoy sequence
            Assert.That(allPsmsArray[5].IsDecoy);

            SpectralLibrarySearchFunction.CalculateSpectralAngles(sl, allPsmsArray, listOfSortedms2Scans, commonParameters);
            Assert.That(allPsmsArray[5].SpectralAngle, Is.EqualTo(0.69).Within(0.01));

        }

        /// <summary>
        /// Tests that spectral library search works correctly with non-specific enzyme search.
        /// This is critical because non-specific search stores PSMs in a different data structure
        /// (fileSpecificPsmsSeparatedByFdrCategory) and spectral angles must be calculated before
        /// FDR analysis since they are used in PEP calculation.
        /// </summary>
        [Test]
        public static void SpectralLibrarySearchWithNonSpecificSearchTest()
        {
            var testDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch");
            var outputDir = Path.Combine(testDir, @"SpectralLibraryNonSpecificSearchTest");

            string library1 = Path.Combine(testDir, @"P16858_target.msp");
            string library2 = Path.Combine(testDir, @"P16858_decoy.msp");
            string fastaDb = Path.Combine(testDir, @"P16858.fasta");
            string spectraFile = Path.Combine(testDir, @"slicedMouse.raw");

            Directory.CreateDirectory(outputDir);

            // Configure search task for non-specific search with spectral library
            var searchTask = new SearchTask();
            searchTask.SearchParameters.SearchType = SearchType.NonSpecific;
            
            // Configure digestion parameters for non-specific search

            DigestionParams d = new DigestionParams(protease: "non-specific", maxMissedCleavages:20, minPeptideLength:5, searchModeType: Omics.Digestion.CleavageSpecificity.None);

            searchTask.CommonParameters = new CommonParameters(
                digestionParams: d);
            searchTask.SearchParameters = new SearchParameters
            {
                SearchType = SearchType.NonSpecific,
            };

            searchTask.RunTask(outputDir,
                new List<DbForTask>
                {
                    new DbForTask(library1, false),
                    new DbForTask(library2, false),
                    new DbForTask(fastaDb, false)
                },
                new List<string> { spectraFile },
                "");

            var results = File.ReadAllLines(Path.Combine(outputDir, @"AllPSMs.psmtsv"));
            var split = results[0].Split('\t');
            
            // Verify spectral angle column exists
            int spectralAngleIndex = Array.IndexOf(split, "Normalized Spectral Angle");
            Assert.That(spectralAngleIndex >= 0, "Normalized Spectral Angle column should exist in output");

            // Verify at least one PSM has a spectral angle calculated
            bool hasSpectralAngle = false;
            for (int i = 1; i < results.Length; i++)
            {
                var columns = results[i].Split('\t');
                if (double.TryParse(columns[spectralAngleIndex], out double spectralAngle) && spectralAngle > 0)
                {
                    hasSpectralAngle = true;
                    break;
                }
            }
            Assert.That(hasSpectralAngle, Is.True, "At least one PSM should have a spectral angle calculated for non-specific search");

            Directory.Delete(outputDir, true);
        }

    }
    
}
