using EngineLayer;
using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using TaskLayer;
using System;
using System.Reflection;
using EngineLayer.DatabaseLoading;
using System.Diagnostics;
using MzLibUtil;
using Omics;
using Omics.Modifications;
using Proteomics;
using Readers;
using Transcriptomics;


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
        [TestCase(SearchType.Classic)]
        [TestCase(SearchType.Modern)]
        public static void TestToleranceExpansion(SearchType searchType)
        {
            // capture warnings
            var originalOut = Console.Out; 
            var originalErr = Console.Error; 
            var sw = new StringWriter(); 
            var listener = new TextWriterTraceListener(sw) { TraceOutputOptions = TraceOptions.None }; 
            Trace.Listeners.Add(listener); 
            Console.SetOut(sw); 
            Console.SetError(sw);

            try
            {
                // set up directories
                string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"ExperimentalDesignCalibrationTest");
                string outputFolder = Path.Combine(unitTestFolder, @"TaskOutput");
                Directory.CreateDirectory(unitTestFolder);
                Directory.CreateDirectory(outputFolder);

                // set up original spectra file (input to calibration)
                string nonCalibratedFilePath = Path.Combine(unitTestFolder, "filename1.mzML");
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML"), nonCalibratedFilePath, true);

                // protein db
                string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

                // run calibration
                CalibrationTask calibrationTask = new();
                calibrationTask.CalibrationParameters.SearchType = searchType;
                calibrationTask.CommonParameters.PrecursorMassTolerance = new PpmTolerance(2);
                calibrationTask.CommonParameters.ProductMassTolerance = new PpmTolerance(2);
                calibrationTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { nonCalibratedFilePath }, "test");


                string expectedCalibratedFileName = Path.GetFileNameWithoutExtension(nonCalibratedFilePath) + "-calib.mzML";
                var expectedCalibratedFilePath = Path.Combine(outputFolder, expectedCalibratedFileName);

                // test file-specific toml written by calibration w/ suggested ppm tolerances
                string expectedTomlName = Path.GetFileNameWithoutExtension(nonCalibratedFilePath) + "-calib.toml";

                Assert.That(File.Exists(Path.Combine(outputFolder, expectedTomlName)));

                var lines = File.ReadAllLines(Path.Combine(outputFolder, expectedTomlName));
                var tolerance = Regex.Match(lines[0], @"\d+\.\d*").Value;
                var tolerance1 = Regex.Match(lines[1], @"\d+\.\d*").Value;

                Assert.That(double.TryParse(tolerance, out double tol));
                Assert.That(double.TryParse(tolerance1, out double tol1));
                Assert.That(lines[0].Contains("PrecursorMassTolerance"));
                Assert.That(lines[1].Contains("ProductMassTolerance"));

                // check that calibrated .mzML exists
                Assert.That(File.Exists(Path.Combine(outputFolder, expectedCalibratedFilePath)));

                // assert warning for wider tolerance occurred
                var output = sw.ToString();
                Assert.That(output, Does.Contain("Could not find enough PSMs to calibrate with; opening up tolerances to"));
                Assert.That(output, Does.Contain("ppm precursor and"));
                Assert.That(output, Does.Contain("ppm product"));

                // clean up
                Directory.Delete(unitTestFolder, true);

            }
            finally
            {
                Trace.Listeners.Remove(listener);
                listener.Flush();
                listener.Close();
                Console.SetOut(originalOut);
                Console.SetError(originalErr);
            }
        }

        [Test]
        [TestCase(SearchType.Classic)]
        [TestCase(SearchType.Modern)]
        public static void CalibrationTestNoPsms(SearchType searchType)
        {
            // set up directories
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"ExperimentalDesignCalibrationTest");
            string outputFolder = Path.Combine(unitTestFolder, @"TaskOutput");
            Directory.CreateDirectory(unitTestFolder);
            Directory.CreateDirectory(outputFolder);

            // set up original spectra file (input to calibration)
            string nonCalibratedFilePath = Path.Combine(unitTestFolder, "filename1.mzML");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML"), nonCalibratedFilePath, true);

            // protein db for a non-matching organism
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\gapdh.fa");

            // set up original experimental design (input to calibration)
            SpectraFileInfo fileInfo = new(nonCalibratedFilePath, "condition", 0, 0, 0);
            _ = ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo> { fileInfo });

            // run calibration
            CalibrationTask calibrationTask = new();
            calibrationTask.CalibrationParameters.SearchType = searchType;
            calibrationTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { nonCalibratedFilePath }, "test");

            // test new experimental design written by calibration
            var newExpDesignPath = Path.Combine(outputFolder, @"ExperimentalDesign.tsv");
            string expectedCalibratedFileName = Path.GetFileNameWithoutExtension(nonCalibratedFilePath) + "-calib.mzML";
            var expectedCalibratedFilePath = Path.Combine(outputFolder, expectedCalibratedFileName);
            var newExperDesign = ExperimentalDesign.ReadExperimentalDesign(newExpDesignPath, new List<string> { expectedCalibratedFilePath }, out var errors);

            Assert.That(errors.Any());

            // clean up
            Directory.Delete(unitTestFolder, true);
        }

        [Test]
        [NonParallelizable]
        public static void CalibrationTooFewMS1DataPoints()
        {
            // set up directories
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"ExperimentalDesignCalibrationTest");
            string outputFolder = Path.Combine(unitTestFolder, @"TaskOutput");
            Directory.CreateDirectory(unitTestFolder);
            Directory.CreateDirectory(outputFolder);

            // set up original spectra file (input to calibration)
            string nonCalibratedFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\mouseOne.mzML");

            // set up original experimental design (input to calibration)
            SpectraFileInfo fileInfo = new(nonCalibratedFilePath, "condition", 0, 0, 0);
            _ = ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo> { fileInfo });

            // protein db for a non-matching organism
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\mouseOne.xml");

            CalibrationTask calibrationTask = new();

            calibrationTask.CommonParameters = new CommonParameters(
                               trimMsMsPeaks: true,
                                              doPrecursorDeconvolution: true);
            var wasCalled = false;

            EventHandler<StringEventArgs> handler = (o, e) => CalibrationWarnHandler(o, e, ref wasCalled);
            MetaMorpheusTask.WarnHandler += handler;

            MetaMorpheusTask.WarnHandler -= (o, e) =>
            {
                wasCalled = true;
                Assert.That(e.S, Does.Contain("Calibration failure! Could not find enough MS1 datapoints."));
            };

            // clean up
            Directory.Delete(unitTestFolder, true);
            MetaMorpheusTask.WarnHandler -= handler;
        }

        private static void CalibrationWarnHandler(object sender, StringEventArgs e, ref bool wasCalled)
        {
            wasCalled = true;
            Assert.That(e.S, Does.Contain("Calibration failure! Could not find enough MS1 datapoints."));
        }

        [Test]
        [TestCase(SearchType.Classic)]
        [TestCase(SearchType.Modern)]
        public static void CalibrationTestLowRes(SearchType searchType)
        {
            CalibrationTask calibrationTask = new CalibrationTask();
            calibrationTask.CalibrationParameters.SearchType = searchType;

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestCalibrationLow");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.fasta");
            Directory.CreateDirectory(outputFolder);

            calibrationTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");
            Assert.That(File.Exists(Path.Combine(outputFolder, @"TaGe_SA_A549_3_snip-calib.mzML")));
            Assert.That(File.Exists(Path.Combine(outputFolder, @"TaGe_SA_A549_3_snip-calib.toml")));
            var lines = File.ReadAllLines(Path.Combine(outputFolder, @"TaGe_SA_A549_3_snip-calib.toml"));
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
        public static void CalibrationTestYeastLowRes()
        {
            CalibrationTask calibrationTask = new CalibrationTask();

            CommonParameters CommonParameters = new(dissociationType: DissociationType.LowCID,
                scoreCutoff: 1);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestCalibrationLow");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            Directory.CreateDirectory(outputFolder);

            calibrationTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");
            Assert.That(File.Exists(Path.Combine(outputFolder, @"SmallCalibratible_Yeast-calib.mzML")));
            Assert.That(File.Exists(Path.Combine(outputFolder, @"SmallCalibratible_Yeast-calib.toml")));
            var lines = File.ReadAllLines(Path.Combine(outputFolder, @"SmallCalibratible_Yeast-calib.toml"));
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
        public static void ExperimentalDesignCalibrationAndSearchWithOneCalibratibleAndOneNoncalibratible()
        {
            // set up output directories
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"ExperimentalDesignCalibrationAndSearch");
            string outputFolder = Path.Combine(unitTestFolder, @"TaskOutput");
            Directory.CreateDirectory(unitTestFolder);
            Directory.CreateDirectory(outputFolder);

            // set up original spectra file (input to calibration)
            string nonCalibratedFilePathOne = Path.Combine(unitTestFolder, "filename1.mzML");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML"), nonCalibratedFilePathOne, true);
            string nonCalibratedFilePathTwo = Path.Combine(unitTestFolder, "filename2.mzML");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.mzML"), nonCalibratedFilePathTwo, true);

            // set up original experimental design (input to calibration)
            SpectraFileInfo fileInfoOne = new(nonCalibratedFilePathOne, "condition1", 0, 0, 0);
            SpectraFileInfo fileInfoTwo = new(nonCalibratedFilePathTwo, "condition2", 0, 0, 0);
            _ = ExperimentalDesign.WriteExperimentalDesignToFile(new List<SpectraFileInfo> { fileInfoOne, fileInfoTwo });

            // set up tasks (calibration + search)
            CalibrationTask calibrationTask = new CalibrationTask();
            SearchTask searchTask = new SearchTask();

            // protein db
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            // run the tasks
            EverythingRunnerEngine a = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("", calibrationTask), ("", searchTask) },
                new List<string> { nonCalibratedFilePathOne, nonCalibratedFilePathTwo },
                new List<DbForTask> { new DbForTask(myDatabase, false) },
                outputFolder);

            a.Run();

            // test to see if quantification ran correctly
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
            Assert.That(errors[0].ToString().Contains(expectedError));
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
            Assert.That(wasCalled);
            

            // clean up
            Directory.Delete(unitTestFolder, true);
        }

        /// <summary>
        /// Tests that calibration using Modern Search produces the expected output files.
        /// This is critical because Modern Search uses an indexed database approach (IndexingEngine + ModernSearchEngine)
        /// rather than the ClassicSearchEngine. We need to verify that the calibration workflow works correctly
        /// when using this alternative search strategy.
        /// </summary>
        [Test]
        public static void CalibrationWithModernSearchTest()
        {
            // set up directories
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"CalibrationModernSearchTest");
            string outputFolder = Path.Combine(unitTestFolder, @"TaskOutput");
            Directory.CreateDirectory(unitTestFolder);
            Directory.CreateDirectory(outputFolder);

            // set up original spectra file (input to calibration)
            string nonCalibratedFilePath = Path.Combine(unitTestFolder, "testfile.mzML");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML"), nonCalibratedFilePath, true);

            // protein db
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            // run calibration with Modern Search
            CalibrationTask calibrationTask = new();
            calibrationTask.CalibrationParameters.SearchType = SearchType.Modern;
            calibrationTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { nonCalibratedFilePath }, "test");

            // test file-specific toml written by calibration w/ suggested ppm tolerances
            string expectedTomlName = Path.GetFileNameWithoutExtension(nonCalibratedFilePath) + "-calib.toml";
            string expectedCalibratedFileName = Path.GetFileNameWithoutExtension(nonCalibratedFilePath) + "-calib.mzML";
            var expectedCalibratedFilePath = Path.Combine(outputFolder, expectedCalibratedFileName);

            Assert.That(File.Exists(Path.Combine(outputFolder, expectedTomlName)), "Calibration toml file should exist");
            Assert.That(File.Exists(expectedCalibratedFilePath), "Calibrated mzML file should exist");

            var lines = File.ReadAllLines(Path.Combine(outputFolder, expectedTomlName));
            var tolerance = Regex.Match(lines[0], @"\d+\.\d*").Value;
            var tolerance1 = Regex.Match(lines[1], @"\d+\.\d*").Value;

            Assert.That(double.TryParse(tolerance, out double tol), "Precursor tolerance should be parseable");
            Assert.That(double.TryParse(tolerance1, out double tol1), "Product tolerance should be parseable");
            Assert.That(lines[0].Contains("PrecursorMassTolerance"), "First line should contain PrecursorMassTolerance");
            Assert.That(lines[1].Contains("ProductMassTolerance"), "Second line should contain ProductMassTolerance");

            // clean up
            Directory.Delete(unitTestFolder, true);
        }

        /// <summary>
        /// Compares calibration results between Classic Search and Modern Search to ensure they produce
        /// comparable results. This is critical to verify that the Modern Search implementation doesn't
        /// introduce significant differences in calibration quality. Both search types should find
        /// similar numbers of PSMs and produce similar mass tolerance recommendations.
        /// </summary>
        [Test]
        public static void CalibrationClassicVsModernSearchComparison()
        {
            // set up directories
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"CalibrationClassicVsModernTest");
            string classicOutputFolder = Path.Combine(unitTestFolder, @"ClassicOutput");
            string modernOutputFolder = Path.Combine(unitTestFolder, @"ModernOutput");
            Directory.CreateDirectory(unitTestFolder);
            Directory.CreateDirectory(classicOutputFolder);
            Directory.CreateDirectory(modernOutputFolder);

            // set up original spectra files (need separate copies to avoid file locking issues)
            string classicFilePath = Path.Combine(unitTestFolder, "classic_test.mzML");
            string modernFilePath = Path.Combine(unitTestFolder, "modern_test.mzML");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML"), classicFilePath, true);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML"), modernFilePath, true);

            // protein db
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            // run calibration with Classic Search
            CalibrationTask classicCalibrationTask = new();
            classicCalibrationTask.CalibrationParameters.SearchType = SearchType.Classic;
            classicCalibrationTask.RunTask(classicOutputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { classicFilePath }, "classic");

            // run calibration with Modern Search
            CalibrationTask modernCalibrationTask = new();
            modernCalibrationTask.CalibrationParameters.SearchType = SearchType.Modern;
            modernCalibrationTask.RunTask(modernOutputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { modernFilePath }, "modern");

            // verify both produced output files
            string classicTomlPath = Path.Combine(classicOutputFolder, "classic_test-calib.toml");
            string modernTomlPath = Path.Combine(modernOutputFolder, "modern_test-calib.toml");

            Assert.That(File.Exists(classicTomlPath), "Classic calibration should produce toml file");
            Assert.That(File.Exists(modernTomlPath), "Modern calibration should produce toml file");

            // read tolerances from both
            var classicLines = File.ReadAllLines(classicTomlPath);
            var modernLines = File.ReadAllLines(modernTomlPath);

            Assert.That(double.TryParse(Regex.Match(classicLines[0], @"\d+\.\d*").Value, out double classicPrecursorTol),
                "Classic precursor tolerance should be parseable");
            Assert.That(double.TryParse(Regex.Match(classicLines[1], @"\d+\.\d*").Value, out double classicProductTol),
                "Classic product tolerance should be parseable");
            Assert.That(double.TryParse(Regex.Match(modernLines[0], @"\d+\.\d*").Value, out double modernPrecursorTol),
                "Modern precursor tolerance should be parseable");
            Assert.That(double.TryParse(Regex.Match(modernLines[1], @"\d+\.\d*").Value, out double modernProductTol),
                "Modern product tolerance should be parseable");

            // tolerances should be very close — both engines search the same data
            // and feed into the same calibration math
            Assert.That(modernPrecursorTol, Is.InRange(classicPrecursorTol * 0.9, classicPrecursorTol * 1.1),
                $"Modern precursor tolerance ({modernPrecursorTol}) should be within 10% of Classic ({classicPrecursorTol})");
            Assert.That(modernProductTol, Is.InRange(classicProductTol * 0.9, classicProductTol * 1.1),
                $"Modern product tolerance ({modernProductTol}) should be within 10% of Classic ({classicProductTol})");

            // clean up
            Directory.Delete(unitTestFolder, true);
        }

        /// <summary>
        /// Verifies that Modern Search calibration throws an informative MetaMorpheusException
        /// when the loaded database contains only non-Protein biopolymers (e.g., RNA).
        /// This covers the proteinList.Count == 0 guard in GenerateIndexes, ensuring users
        /// get a clear error message directing them to use Classic Search for non-protein analytes.
        /// </summary>
        [Test]
        public static void ModernSearchCalibration_AllNonProteinBiopolymers_ThrowsInformativeException()
        {
            CalibrationTask calibrationTask = new();
            calibrationTask.CalibrationParameters.SearchType = SearchType.Modern;

            // Set up the private fields that GenerateIndexes depends on via reflection.
            // _proteinList contains only RNA (non-Protein IBioPolymer) — this triggers the exception.
            var rnaList = new List<IBioPolymer>
            {
                new RNA("GUACUG", "rna_1"),
                new RNA("AACUGCUAG", "rna_2")
            };

            SetPrivateField(calibrationTask, "_proteinList", rnaList);
            SetPrivateField(calibrationTask, "_taskId", "test");
            SetPrivateField(calibrationTask, "_variableModifications", new List<Modification>());
            SetPrivateField(calibrationTask, "_fixedModifications", new List<Modification>());
            SetPrivateField(calibrationTask, "_dbFilenameList", new List<DbForTask>());

            // Invoke the private GenerateIndexes method directly
            MethodInfo generateIndexes = typeof(CalibrationTask)
                .GetMethod("GenerateIndexes", BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.That(generateIndexes, Is.Not.Null, "GenerateIndexes method should exist");

            CommonParameters combinedParams = new();

            // The reflection call wraps the inner exception in TargetInvocationException
            var ex = Assert.Throws<TargetInvocationException>(() =>
                generateIndexes.Invoke(calibrationTask, new object[] { combinedParams, "testfile" }));

            // Verify it's the expected MetaMorpheusException with a helpful message
            Assert.That(ex.InnerException, Is.TypeOf<MetaMorpheusException>());
            Assert.That(ex.InnerException.Message, Does.Contain("Modern Search calibration requires protein sequences"));
            Assert.That(ex.InnerException.Message, Does.Contain("2 non-protein biopolymer(s) were filtered out"),
                "Exception should report the count of filtered-out biopolymers");
            Assert.That(ex.InnerException.Message, Does.Contain("Use Classic Search calibration"),
                "Exception should suggest the alternative search type");
        }

        /// <summary>
        /// Verifies that Modern Search calibration emits a warning when the database contains
        /// a mix of Protein and non-Protein biopolymers (e.g., Protein + RNA). The non-protein
        /// entries are silently excluded from the search index, and the user should be warned.
        /// This covers the proteinList.Count &lt; _proteinList.Count branch in GenerateIndexes.
        /// </summary>
        [Test]
        public static void ModernSearchCalibration_MixedBiopolymers_WarnsAboutExcludedEntries()
        {
            CalibrationTask calibrationTask = new();
            calibrationTask.CalibrationParameters.SearchType = SearchType.Modern;

            // Create a mixed list: one real Protein (needed for IndexingEngine to succeed) and one RNA.
            // The RNA entry will be filtered out, triggering the warning.
            var protein = new Protein("MKWVTFISLLLLFSSAYSR", "P00001");
            var rna = new RNA("GUACUG", "rna_1");
            var mixedList = new List<IBioPolymer> { protein, rna };

            SetPrivateField(calibrationTask, "_proteinList", mixedList);
            SetPrivateField(calibrationTask, "_taskId", "test");
            SetPrivateField(calibrationTask, "_variableModifications", new List<Modification>());
            SetPrivateField(calibrationTask, "_fixedModifications", new List<Modification>());
            SetPrivateField(calibrationTask, "_dbFilenameList", new List<DbForTask>());

            // Subscribe to WarnHandler to capture the warning message
            string capturedWarning = null;
            EventHandler<StringEventArgs> handler = (_, e) => capturedWarning = e.S;
            MetaMorpheusTask.WarnHandler += handler;

            try
            {
                // Call GenerateIndexes — the warning fires before IndexingEngine runs.
                // IndexingEngine may throw due to minimal setup, but the warning should already be captured.
                MethodInfo generateIndexes = typeof(CalibrationTask)
                    .GetMethod("GenerateIndexes", BindingFlags.NonPublic | BindingFlags.Instance);
                Assert.That(generateIndexes, Is.Not.Null, "GenerateIndexes method should exist");

                CommonParameters combinedParams = new();

                try
                {
                    generateIndexes.Invoke(calibrationTask, new object[] { combinedParams, "testfile" });
                }
                catch (TargetInvocationException)
                {
                    // IndexingEngine may fail with minimal test setup — that's fine.
                    // The warning we're testing fires before IndexingEngine runs.
                }

                // Verify the warning was emitted with the correct message
                Assert.That(capturedWarning, Is.Not.Null, "A warning should be emitted when non-protein biopolymers are excluded");
                Assert.That(capturedWarning, Does.Contain("1 non-protein biopolymer(s) were excluded from the search index"),
                    "Warning should report the exact count of excluded entries");
                Assert.That(capturedWarning, Does.Contain("Only protein sequences are supported by Modern Search"),
                    "Warning should explain the limitation");
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
            }
        }

        /// <summary>
        /// Helper to set private/internal fields on CalibrationTask via reflection.
        /// Used by tests that need to control internal state for targeted unit testing.
        /// </summary>
        private static void SetPrivateField(object obj, string fieldName, object value)
        {
            FieldInfo field = obj.GetType().GetField(fieldName,
                BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.That(field, Is.Not.Null, $"Field '{fieldName}' should exist on {obj.GetType().Name}");
            field.SetValue(obj, value);
        }

    }
}