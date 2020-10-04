using EngineLayer;
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
        public static void TestScanStreaming()
        {
            string path = @"C:\Data\Yeast\09-04-18_EcoliSpikeInSingleShot1x.mzML";

            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            ScanStreamer s = new ScanStreamer(path, new CommonParameters());

            List<string> output = new List<string>();

            Parallel.ForEach(Partitioner.Create(0, 1000000), new ParallelOptions { MaxDegreeOfParallelism = -1 },
                (range, loopState) =>
                {
                    Ms2ScanWithSpecificMass scan;

                    while ((scan = s.NextScan()) != null)
                    {
                        lock (output)
                        {
                            output.Add(scan.OneBasedScanNumber.ToString());
                        }
                    }
                });

            //Ms2ScanWithSpecificMass scan;
            //while ((scan = s.NextScan()) != null)
            //{
            //    scan = s.NextScan();
            //}

            stopwatch.Stop();
            stopwatch.Restart();

            var ff = Mzml.LoadAllStaticData(path);
            var scans = MetaMorpheusTask.GetMs2Scans(ff, path, new CommonParameters()).ToList();

            foreach (var scan in scans)
            {

            }

            stopwatch.Stop();
        }
    }
}