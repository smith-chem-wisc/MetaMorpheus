using System;
using System.Collections.Generic;
using EngineLayer.ClassicSearch;
using EngineLayer;
using System.Reflection;
using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using EngineLayer.Calibration;
using TaskLayer;
using EngineLayer.FdrAnalysis;
using MzLibUtil;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class CalibrationRetainsProperTolerances
    {
        [Test]
        public void EnsureTolerancesAreCorrectThroughout()
        {
            // set up directories and input data
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"CaliConsistencyTest");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);
            Directory.CreateDirectory(outputFolder);

            CalibrationTask calibrationTask = new();
            string nonCalibratedFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");

            // define expected values
            int expectedCalibrationCount = 2;
            int expectedClassicSearchCount = 3;
            int expectedDataPointAcquisitionCount = 3;
            int expectedFdrAnalysisEngineCount = 3;

            double[] expectedPrecursorTolerances = [10, 8.6, 7.3];
            double[] expectedProductTolerances = [30, 14.4, 13.9];

            int classicSearchEngineCount = 0;
            int calibrationEngineCount = 0;
            int dataPointAcquisitionCount = 0;
            int fdrAnalysisEngineCount = 0;

            List<(double precursor, double product, double expectedPrecursor, double expectedProduct)> classicSearchResults = new();
            List<(double precursor, double product, double expectedPrecursor, double expectedProduct)> calibrationResults = new();
            List<(double precursor, double product, double expectedPrecursor, double expectedProduct)> dataPointAcquisitionResults = new();

            // In one Calibration task, the ClassicSearchEngine, DataPointAcquisitionEngine, and FdrAnalysisEngine all run three times,
            // Each time, calibration is refined and tolerances get progressively smaller. 
            // This event handler listens for engine calls, and records what tolerances are used on each call.
            // This is compared to the expected values, which start with the initial default tolerances (10 precursor, 30 product),
            // and then get smaller with each calibration round in a way which is dependent on the data being calibrated
            EventHandler<SingleEngineFinishedEventArgs> finishedEngineEventHandler = (sender, e) =>
            {
                switch (sender)
                {
                    case null:
                        Assert.Fail("Cannot check engine type");
                        break;
                    case ClassicSearchEngine cse:
                        {
                            classicSearchEngineCount++;

                            // Ensure precursor tol and mass diff acceptor agree
                            var precursor = cse.CommonParameters.PrecursorMassTolerance.Value;
                            var searchMode = GetSearchMode(cse);
                            var massDiffAcceptorString = searchMode.ToProseString();
                            var precursorFromAcceptor = double.Parse(massDiffAcceptorString.Split(' ')[0]);
                            Assert.That(precursor, Is.EqualTo(precursorFromAcceptor));

                            var product = cse.CommonParameters.ProductMassTolerance.Value;

                            var expectedPrecursor = expectedPrecursorTolerances[classicSearchEngineCount - 1];
                            var expectedProduct = expectedProductTolerances[classicSearchEngineCount - 1];

                            classicSearchResults.Add((precursor, product, expectedPrecursor, expectedProduct));
                        }
                        break;
                    // CalibrationEngine only twice, once after the initial round of calibration, once after the second round of calibration
                    case CalibrationEngine cale:
                        {
                            calibrationEngineCount++;

                            var precursor = cale.CommonParameters.PrecursorMassTolerance.Value;
                            var product = cale.CommonParameters.ProductMassTolerance.Value;

                            var expectedPrecursor = expectedPrecursorTolerances[calibrationEngineCount];
                            var expectedProduct = expectedProductTolerances[calibrationEngineCount];

                            calibrationResults.Add((precursor, product, expectedPrecursor, expectedProduct));
                        }
                        break;
                    case DataPointAcquisitionEngine dpae:
                        {
                            dataPointAcquisitionCount++;
                            var precursor = ((Tolerance)dpae.GetType().GetField("PrecursorMassTolerance", BindingFlags.NonPublic | BindingFlags.Instance)!
                                .GetValue(dpae)!).Value;
                            var product = ((Tolerance)dpae.GetType().GetField("ProductMassTolerance", BindingFlags.NonPublic | BindingFlags.Instance)!
                                    .GetValue(dpae)!).Value;

                            var expectedPrecursor = expectedPrecursorTolerances[dataPointAcquisitionCount - 1];
                            var expectedProduct = expectedProductTolerances[dataPointAcquisitionCount - 1];

                            dataPointAcquisitionResults.Add((precursor, product, expectedPrecursor, expectedProduct));
                        }
                        break;
                    case FdrAnalysisEngine fdre:
                        fdrAnalysisEngineCount++;
                        break;

                    default:
                        Assert.Fail("Unexpected Engine Type");
                        break;
                }
            };

            // Ensure we catch the end of all engines running and record the tolerances
            MetaMorpheusEngine.FinishedSingleEngineHandler += finishedEngineEventHandler;

            // run calibration
            calibrationTask.RunTask(outputFolder, [new DbForTask(myDatabase, false)], [nonCalibratedFilePath], "test");

            // Unsubscribe from the event
            MetaMorpheusEngine.FinishedSingleEngineHandler -= finishedEngineEventHandler;

            Assert.That(calibrationEngineCount, Is.EqualTo(expectedCalibrationCount));
            Assert.That(classicSearchEngineCount, Is.EqualTo(expectedClassicSearchCount));
            Assert.That(dataPointAcquisitionCount, Is.EqualTo(expectedDataPointAcquisitionCount));
            Assert.That(fdrAnalysisEngineCount, Is.EqualTo(expectedFdrAnalysisEngineCount));

            foreach (var result in classicSearchResults)
            {
                Assert.That(result.precursor, Is.EqualTo(result.expectedPrecursor));
                Assert.That(result.product, Is.EqualTo(result.expectedProduct));
            }

            foreach (var result in calibrationResults)
            {
                Assert.That(result.precursor, Is.EqualTo(result.expectedPrecursor));
                Assert.That(result.product, Is.EqualTo(result.expectedProduct));
            }

            foreach (var result in dataPointAcquisitionResults)
            {
                Assert.That(result.precursor, Is.EqualTo(result.expectedPrecursor));
                Assert.That(result.product, Is.EqualTo(result.expectedProduct));
            }

            Assert.That(calibrationEngineCount, Is.EqualTo(expectedCalibrationCount));
            Assert.That(classicSearchEngineCount, Is.EqualTo(expectedClassicSearchCount));
            Assert.That(dataPointAcquisitionCount, Is.EqualTo(expectedDataPointAcquisitionCount));
            Assert.That(fdrAnalysisEngineCount, Is.EqualTo(expectedFdrAnalysisEngineCount));

            foreach (var result in classicSearchResults)
            {
                Assert.That(result.precursor, Is.EqualTo(result.expectedPrecursor));
                Assert.That(result.product, Is.EqualTo(result.expectedProduct));
            }

            foreach (var result in calibrationResults)
            {
                Assert.That(result.precursor, Is.EqualTo(result.expectedPrecursor));
                Assert.That(result.product, Is.EqualTo(result.expectedProduct));
            }

            foreach (var result in dataPointAcquisitionResults)
            {
                Assert.That(result.precursor, Is.EqualTo(result.expectedPrecursor));
                Assert.That(result.product, Is.EqualTo(result.expectedProduct));
            }
        }

        private static MassDiffAcceptor GetSearchMode(ClassicSearchEngine engine)
        {
            // Get the type of the ClassicSearchEngine
            Type type = typeof(ClassicSearchEngine);

            // Get the private field 'SearchMode'
            FieldInfo fieldInfo = type.GetField("SearchMode", BindingFlags.NonPublic | BindingFlags.Instance);

            // Get the value of the 'SearchMode' field from the instance of ClassicSearchEngine
            return (MassDiffAcceptor)fieldInfo!.GetValue(engine);
        }
    }
}
