using EngineLayer;
using IO.Mgf;
using IO.MzML;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class MsDataFileTest
    {
        [OneTimeSetUp]
        public static void Setup()
        {
            Environment.CurrentDirectory = TestContext.CurrentContext.TestDirectory;
        }

        [Test]
        public static void TestLoadAndRunMgf()
        {
            //The purpose of this test is to ensure that mgfs can be run without crashing.
            //Whenever a new feature is added that may require things an mgf does not have,
            //there should be a check that prevents mgfs from using that feature.
            string mgfName = @"TestData\ok.mgf";
            string xmlName = @"TestData\okk.xml";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestLoadAndRunMgf");

            SearchTask task1 = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DoQuantification = true
                }
            };
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)>
            {
                ("task1", task1),
            };
            //run!

            var engine = new EverythingRunnerEngine(taskList, new List<string> { mgfName }, new List<DbForTask> { new DbForTask(xmlName, false) }, outputFolder);
            engine.Run();
            //Just don't crash! There should also be at least one psm at 1% FDR, but can't check for that.
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestCompressionDecompression()
        {
            string testInputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"CompressionTest");
            DirectoryInfo testDirectory = new DirectoryInfo(testInputFolder);
            MyFileManager.CompressDirectory(testDirectory);

            foreach (FileInfo file in testDirectory.GetFiles())
            {
                Assert.AreEqual(".gz", file.Extension);
            }

            MyFileManager.DecompressDirectory(testDirectory);

            foreach (FileInfo file in testDirectory.GetFiles())
            {
                Assert.AreNotEqual(".gz", file.Extension);
            }
        }

        [Test]
        public static void TestMs2ScanWithSpecificMass()
        {
            Ms2ScanWithSpecificMass scanB = new Ms2ScanWithSpecificMass(
                new MsDataScan(
                    new MzSpectrum(new double[] { }, new double[] { }, false),
                    2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null),
                100, 1, null, new CommonParameters(), null);

            var closestExperimentalMassB = scanB.GetClosestExperimentalIsotopicEnvelope(10);

            Assert.IsNull(closestExperimentalMassB);
        }

        [Test]
        [TestCase(@"TestData\SmallCalibratible_Yeast.mzML")]
        public static void TestScanStreaming(string file)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, file);

            // need to use only 1 thread in this unit test because the order of Ms2ScanWithSpecificMass objects 
            // in the foreach loop below matters. the MetaMorpheusTask.GetMs2Scans method is parallelized
            // and does not return scans in a particular order
            var commonParam = new CommonParameters(maxThreadsToUsePerFile: 1);

            ScanStreamer s = new ScanStreamer(path, commonParam);
            var filter = MyFileManager.GetFilterParamsFromCommonParams(commonParam);

            MsDataFile staticFile = null;

            var ext = Path.GetExtension(path).ToLowerInvariant();
            if (ext == ".mzml")
            {
                staticFile = Mzml.LoadAllStaticData(path, filter);
            }
            else if (ext == ".raw")
            {
                staticFile = ThermoRawFileReader.LoadAllStaticData(path, filter);
            }
            else if (ext == ".mgf")
            {
                staticFile = Mgf.LoadAllStaticData(path, filter);
            }

            var staticScans = MetaMorpheusTask.GetMs2Scans(staticFile, path, new CommonParameters()).ToList();

            foreach (Ms2ScanWithSpecificMass staticScan in staticScans)
            {
                Ms2ScanWithSpecificMass streamedScan = s.NextScan();

                Assert.That(staticScan.PrecursorMass == streamedScan.PrecursorMass);
                Assert.That(staticScan.PrecursorCharge == streamedScan.PrecursorCharge);
                Assert.That(staticScan.OneBasedPrecursorScanNumber == streamedScan.OneBasedPrecursorScanNumber);
                Assert.That(staticScan.OneBasedScanNumber == streamedScan.OneBasedScanNumber);

                Assert.That(staticScan.ExperimentalFragments.Length == streamedScan.ExperimentalFragments.Length);

                for (int i = 0; i < staticScan.ExperimentalFragments.Length; i++)
                {
                    var staticFragment = staticScan.ExperimentalFragments[i];
                    var dynamicFragment = streamedScan.ExperimentalFragments[i];

                    Assert.That(staticFragment.MonoisotopicMass == dynamicFragment.MonoisotopicMass);
                    Assert.That(staticFragment.Charge == dynamicFragment.Charge);
                }
            }

            Assert.That(s.NextScan() == null);


            //string path = @"C:\Data\Yeast\09-04-18_EcoliSpikeInSingleShot1x.mzML";

            //Stopwatch stopwatch = new Stopwatch();
            //stopwatch.Start();

            //List<string> output = new List<string>();

            //Parallel.ForEach(Partitioner.Create(0, 1000000), new ParallelOptions { MaxDegreeOfParallelism = -1 },
            //    (range, loopState) =>
            //    {
            //        Ms2ScanWithSpecificMass scan;

            //        while ((scan = s.NextScan()) != null)
            //        {
            //            lock (output)
            //            {
            //                output.Add(scan.OneBasedScanNumber.ToString());
            //            }
            //        }
            //    });

            //stopwatch.Stop();
        }
    }
}