using Easy.Common.Extensions;
using EngineLayer;
using EngineLayer.Calibration;
using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text.RegularExpressions;
using TaskLayer;
using ThermoFisher.CommonCore.Data;
using ThermoFisher.CommonCore.Data.Business;
using static iText.StyledXmlParser.Jsoup.Select.Evaluator;

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
            SpectraFileInfo fileInfo = new(nonCalibratedFilePath, "condition", 0, 0, 0);
            _ = ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo> { fileInfo });

            // run calibration
            CalibrationTask calibrationTask = new();
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

            CommonParameters CommonParameters = new(dissociationType: DissociationType.LowCID,
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
            SpectraFileInfo fileInfo = new(nonCalibratedFilePath, "condition", 0, 0, 0);
            _ = ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo> { fileInfo });

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

        [Test]
        [TestCase("ExpDesFileNotFound","small.mzML", "Experimental design file not found!")]
        [TestCase("WrongNumberOfCells", "small.mzML", "Error: The experimental design was not formatted correctly. Expected 5 cells, but found 4 on line 2")]
        [TestCase("BioRepNotInteger", "small.mzML", "Error: The experimental design was not formatted correctly. The biorep on line 2 is not an integer")]
        [TestCase("FractionNotInteger", "small.mzML", "Error: The experimental design was not formatted correctly. The fraction on line 2 is not an integer")]
        [TestCase("TechRepNotInt", "small.mzML", "Error: The experimental design was not formatted correctly. The techrep on line 2 is not an integer")]
        [TestCase("mzMLmissing", "small.mzML", "Error: The experimental design did not contain the file(s):")] //tough to check b/c local path not translated to appveyor
        public static void TestExperimentalDesignErrors(string experimentalFolder, string rawFile, string expectedError)
        {
            string experimentalDesignPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData", @"TestExperimentalDesign", experimentalFolder, "ExperimentalDesign.tsv");
            List<string> rawFilePaths = new() { Path.Combine(experimentalDesignPath, rawFile) };
            _ = ExperimentalDesign.ReadExperimentalDesign(experimentalDesignPath, rawFilePaths, out var errors);
            Assert.IsTrue(errors[0].ToString().Contains(expectedError));
        }

        [Test]
        public static void TestWriteNewExperimentalDesignFileDuringCalibration()
        {
            string badExperimentalDesignPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData", @"TestExperimentalDesign", "WrongNumberOfCells", "ExperimentalDesign.tsv");

            // set up directories
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"ExperimentalDesignCalibrationTest");
            string outputFolder = Path.Combine(unitTestFolder, @"TaskOutput");
            Directory.CreateDirectory(unitTestFolder);
            Directory.CreateDirectory(outputFolder);
            
            // set up original spectra file (input to calibration)
            string nonCalibratedFilePath = Path.Combine(unitTestFolder, "filename1.mzML");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML"), nonCalibratedFilePath, true);

            // set up original BAD experimental design (input to calibration)
            string experimentalDesignPath = Path.Combine(unitTestFolder, "ExperimentalDesign.tsv");
            File.Copy(badExperimentalDesignPath, experimentalDesignPath, true);
            
            // protein db
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            // run calibration
            CalibrationTask calibrationTask = new CalibrationTask();


            bool wasCalled = false;
            MetaMorpheusTask.WarnHandler += (o, e) => wasCalled = true;
            
            calibrationTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { nonCalibratedFilePath }, "test");
            
            //The original experimental design file is bad so we expect Warn event in "WriteNewExperimentalDesignFile"
            Assert.IsTrue(wasCalled);
            

            // clean up
            Directory.Delete(unitTestFolder, true);
        }

    }
}