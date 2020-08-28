using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using IO.MzML;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class TestTopDown
    {
        [Test]
        public static void TestClassicSearchEngineTopDown()
        {
            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(protease: "top-down"),
                scoreCutoff: 1,
                assumeOrphanPeaksAreZ1Fragments: false);

            MetaMorpheusTask.DetermineAnalyteType(CommonParameters);

            // test output file name (should be proteoform and not peptide)
            Assert.That(GlobalVariables.AnalyteType == "Proteoform");

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein>
            {
                new Protein("MPKVYSYQEVAEHNGPENFWIIIDDKVYDVSQFKDEHPGGDEIIMDLGGQDATESFVDIGHSDEALRLLKGLYIGDVDKTSERVSVEKVSTSENQSKGSGTLVVILAILMLGVAYYLLNE", "P40312")
            };

            var myMsDataFile = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\slicedTDYeast.mzML"));

            var searchMode = new SinglePpmAroundZeroSearchMode(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null, proteinList, searchMode, CommonParameters, null, new List<string>()).Run();

            var psm = allPsmsArray.Where(p => p != null).FirstOrDefault();
            Assert.That(psm.MatchedFragmentIons.Count == 47);
        }

        [Test]
        public static void TestModernSearchEngineTopDown()
        {
            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(protease: "top-down"),
                scoreCutoff: 1,
                assumeOrphanPeaksAreZ1Fragments: false);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein>
            {
                new Protein("MPKVYSYQEVAEHNGPENFWIIIDDKVYDVSQFKDEHPGGDEIIMDLGGQDATESFVDIGHSDEALRLLKGLYIGDVDKTSERVSVEKVSTSENQSKGSGTLVVILAILMLGVAYYLLNE", "P40312")
            };

            var myMsDataFile = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\slicedTDYeast.mzML"));

            var searchMode = new SinglePpmAroundZeroSearchMode(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse, CommonParameters,
                null, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, CommonParameters, null, searchMode, 0, new List<string>()).Run();

            var psm = allPsmsArray.Where(p => p != null).FirstOrDefault();
            Assert.That(psm.MatchedFragmentIons.Count == 47);
        }

        //[Test]
        //public static void FakeTest()
        //{
        //    string spectraFile = @"C:\Users\rmillikin\Desktop\20180913_JMD_MYOGLOBIN_MS2.RAW";

        //    var file = ThermoRawFileReader.LoadAllStaticData(spectraFile);

        //    var scans = file.GetAllScansList().Where(p => p.OneBasedScanNumber == 1269 || p.OneBasedScanNumber == 1270).ToList();

        //    List<string> output = new List<string>();
        //    foreach (var scan in scans)
        //    {
        //        output.Add("BEGIN IONS");
        //        output.Add("TITLE=TEST");
        //        output.Add("PEPMASS=" + scan.SelectedIonMZ);
        //        output.Add("CHARGE=" + scan.SelectedIonChargeStateGuess + "+");
        //        output.Add("RTINSECONDS=" + (int)(scan.RetentionTime * 60));
        //        output.Add("SCANS=" + scan.OneBasedScanNumber);

        //        for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
        //        {
        //            output.Add(scan.MassSpectrum.XArray[i].ToString("F4") + " " + scan.MassSpectrum.YArray[i].ToString("F1"));
        //        }

        //        output.Add("END IONS");
        //    }

        //    string outputFileName = @"C:\Users\rmillikin\Desktop\converted.mgf";
        //    File.WriteAllLines(outputFileName, output);
        //    var temp = IO.Mgf.Mgf.LoadAllStaticData(outputFileName);
        //}

        [Test]
        public static void TestTruncationSearch()
        {
            var temp = ThermoRawFileReader.LoadAllStaticData(@"C:\Data\TD Yeast\05-26-17_B7A_yeast_td_fract7_rep1.raw");
            var ok = temp.GetAllScansList().Where(p => p.OneBasedScanNumber == 1269 || p.OneBasedScanNumber == 1270).ToList();

            var msDataFile = new MsDataFile(ok.ToArray(), new SourceFile(null, null, null, null, null));

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(
                msDataFile,
                @"C:\Data\TD Yeast\05-26-17_B7A_yeast_td_fract7_rep1.mzml",
                false);

            var mzml = Mzml.LoadAllStaticData(@"C:\Data\TD Yeast\05-26-17_B7A_yeast_td_fract7_rep1.mzml");



            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(protease: "top-down"),
                scoreCutoff: 5,
                assumeOrphanPeaksAreZ1Fragments: false);

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein>
            {
                new Protein("MAKSKNHTAHNQTRKAHRNGIKKPKTYKYPSLKGVDPKFRRNHKHALHGTAKALAAAKK", "P05747")
            };

            var myMsDataFile = Mzml.LoadAllStaticData(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\slicedTDYeastTruncation.mzML"));

            var searchMode = new SinglePpmAroundZeroSearchMode(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            var indexEngine = new IndexingEngine(proteinList, variableModifications, fixedModifications, null, null, null, 1, DecoyType.Reverse, CommonParameters,
                null, 30000, false, new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>());
            var indexResults = (IndexingResults)indexEngine.Run();

            new ModernSearchEngine(allPsmsArray, listOfSortedms2Scans, indexResults.PeptideIndex, indexResults.FragmentIndex, 0, CommonParameters, null, searchMode, 0, new List<string>()).Run();

            var psm = allPsmsArray.Where(p => p != null).FirstOrDefault();
            Assert.That(psm.MatchedFragmentIons.Count == 47);
        }
    }
}
