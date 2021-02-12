using EngineLayer;
using FlashLFQ;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using TaskLayer;

namespace Test
{
    internal class CalibrationTests
    {
        [Test]
        [TestCase("filename1.mzML")]
        [TestCase("filename1.1.mzML")]
        public static void ExperimentalDesignCalibrationTest(string nonCalibratedFile)
        {
            // set up directories
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"ExperimentalDesignCalibrationTest");
            string outputFolder = Path.Combine(unitTestFolder, @"TaskOutput");
            Directory.CreateDirectory(unitTestFolder);
            Directory.CreateDirectory(outputFolder);

            // set up original spectra file (input to calibration)
            string nonCalibratedFilePath = Path.Combine(unitTestFolder, nonCalibratedFile);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML"), nonCalibratedFilePath, true);

            // protein db
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            // set up original experimental design (input to calibration)
            SpectraFileInfo fileInfo = new SpectraFileInfo(nonCalibratedFilePath, "condition", 0, 0, 0);
            var experimentalDesignFilePath = ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo> { fileInfo });

            // run calibration
            CalibrationTask calibrationTask = new CalibrationTask();
            calibrationTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { nonCalibratedFilePath }, "test");

            // test new experimental design written by calibration
            var newExpDesignPath = Path.Combine(outputFolder, @"ExperimentalDesign.tsv");
            string expectedCalibratedFileName = Path.GetFileNameWithoutExtension(nonCalibratedFilePath) + "-calib.mzML";
            var expectedCalibratedFilePath = Path.Combine(outputFolder, expectedCalibratedFileName);
            var newExperDesign = ExperimentalDesign.ReadExperimentalDesign(newExpDesignPath, new List<string> { expectedCalibratedFilePath }, out var errors);
            
            Assert.That(!errors.Any());
            Assert.That(newExperDesign.Count == 1);

            // test file-specific toml written by calibration w/ suggested ppm tolerances
            string expectedTomlName = Path.GetFileNameWithoutExtension(nonCalibratedFilePath) + "-calib.toml";

            Assert.That(File.Exists(Path.Combine(outputFolder, expectedTomlName)));

            var lines = File.ReadAllLines(Path.Combine(outputFolder, expectedTomlName));
            var tolerance = Regex.Match(lines[0], @"\d+\.\d*").Value;
            var tolerance1 = Regex.Match(lines[1], @"\d+\.\d*").Value;

            Assert.That(double.TryParse(tolerance, out double tol) == true);
            Assert.That(double.TryParse(tolerance1, out double tol1) == true);
            Assert.That(lines[0].Contains("PrecursorMassTolerance"));
            Assert.That(lines[1].Contains("ProductMassTolerance"));

            // check that calibrated .mzML exists
            Assert.That(File.Exists(Path.Combine(outputFolder, expectedCalibratedFilePath)));

            // clean up
            Directory.Delete(unitTestFolder, true);
        }

        [Test]
        public static void CalibrationTestLowRes()
        {
            CalibrationTask calibrationTask = new CalibrationTask();

            EngineLayer.CommonParameters CommonParameters = new EngineLayer.CommonParameters
                (dissociationType: MassSpectrometry.DissociationType.LowCID,
                scoreCutoff: 1);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestCalibrationLow");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\sliced_b6.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\LowResSnip_B6_mouse_11700_117500pruned.xml");
            Directory.CreateDirectory(outputFolder);

            calibrationTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");
            Assert.That(File.Exists(Path.Combine(outputFolder, @"sliced_b6-calib.mzML")));
            Assert.That(File.Exists(Path.Combine(outputFolder, @"sliced_b6-calib.toml")));
            var lines = File.ReadAllLines(Path.Combine(outputFolder, @"sliced_b6-calib.toml"));
            var tolerance = Regex.Match(lines[0], @"\d+\.\d*").Value;
            var tolerance1 = Regex.Match(lines[1], @"\d+\.\d*").Value;
            Assert.That(double.TryParse(tolerance, out double tol) == true);
            Assert.That(double.TryParse(tolerance1, out double tol1) == true);
            Assert.That(lines[0].Contains("PrecursorMassTolerance"));
            Assert.That(lines[1].Contains("ProductMassTolerance"));
            Directory.Delete(outputFolder, true);
            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Task Settings"), true);
        }

        [Test]
        [TestCase("filename1.1.mzML")]
        public static void ExperimentalDesignCalibrationAndSearch(string nonCalibratedFile)
        {
            // set up output directories
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"ExperimentalDesignCalibrationAndSearch");
            string outputFolder = Path.Combine(unitTestFolder, @"TaskOutput");
            Directory.CreateDirectory(unitTestFolder);
            Directory.CreateDirectory(outputFolder);

            // set up original spectra file (input to calibration)
            string nonCalibratedFilePath = Path.Combine(unitTestFolder, nonCalibratedFile);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML"), nonCalibratedFilePath, true);

            // set up original experimental design (input to calibration)
            SpectraFileInfo fileInfo = new SpectraFileInfo(nonCalibratedFilePath, "condition", 0, 0, 0);
            var experimentalDesignFilePath = ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo> { fileInfo });

            // set up tasks (calibration + search)
            CalibrationTask calibrationTask = new CalibrationTask();
            SearchTask searchTask = new SearchTask();

            // protein db
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            // run the tasks
            EverythingRunnerEngine a = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("", calibrationTask), ("", searchTask) }, 
                new List<string> { nonCalibratedFilePath }, 
                new List<DbForTask> { new DbForTask(myDatabase, false) },
                outputFolder);

            a.Run();

            // test to see if quantification ran correctly
            Assert.That(File.Exists(Path.Combine(outputFolder, @"AllQuantifiedPeptides.tsv")));
            Assert.That(File.Exists(Path.Combine(outputFolder, @"ExperimentalDesign.tsv")));

            // clean up
            Directory.Delete(unitTestFolder, true);
        }
    }
}