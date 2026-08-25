using EngineLayer;
using EngineLayer.DatabaseLoading;
using GuiFunctions;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Threading;
using System.Windows;
using System.Windows.Controls;
using TaskLayer;

namespace Test.MetaDraw
{
    [TestFixture]
    [Apartment(ApartmentState.STA)]
    public class MetaDrawExportRegressionTests
    {
        [SetUp]
        public void SetUp()
        {
            MetaDrawSettings.ResetSettings();
        }

        [Test]
        public void ExportPlot_WhenStationarySequenceDrawingIsDisabled_DoesNotThrow()
        {
            WithLoadedMetaDraw((logic, psm, outputDirectory) =>
            {
                MetaDrawSettings.DrawStationarySequence = false;
                var plotView = new OxyPlot.Wpf.PlotView { Name = "plotView" };
                var parentChildScanPlotsView = new ParentChildScanPlotsView();
                var stationaryCanvas = new Canvas();

                logic.DisplaySequences(stationaryCanvas, null, null, psm);
                logic.DisplaySpectrumMatch(plotView, psm, parentChildScanPlotsView, out var displayErrors);

                Assert.That(logic.StationarySequence, Is.Null);
                Assert.That(displayErrors, Is.Null.Or.Empty);
                Assert.DoesNotThrow(() => logic.ExportPlot(
                    plotView,
                    stationaryCanvas,
                    new List<SpectrumMatchFromTsv> { psm },
                    parentChildScanPlotsView,
                    outputDirectory,
                    out _));
            });
        }

        [Test]
        public void ExportPlot_WhenPsmIsAmbiguous_DoesNotThrow()
        {
            WithLoadedMetaDraw((logic, psm, outputDirectory) =>
            {
                var ambiguousPsm = psm.ReplaceFullSequence(psm.FullSequence + "|" + psm.FullSequence);
                var plotView = new OxyPlot.Wpf.PlotView { Name = "plotView" };
                var parentChildScanPlotsView = new ParentChildScanPlotsView();
                var stationaryCanvas = new Canvas();

                logic.DisplaySequences(stationaryCanvas, null, null, ambiguousPsm);
                logic.DisplaySpectrumMatch(plotView, ambiguousPsm, parentChildScanPlotsView, out var displayErrors);

                Assert.That(logic.StationarySequence, Is.Null);
                Assert.That(displayErrors, Is.Null.Or.Empty);
                Assert.DoesNotThrow(() => logic.ExportPlot(
                    plotView,
                    stationaryCanvas,
                    new List<SpectrumMatchFromTsv> { ambiguousPsm },
                    parentChildScanPlotsView,
                    outputDirectory,
                    out _));
            });
        }

        [Test]
        public void ExportSequenceCoverage_WhenCanvasesHaveNoSize_DoesNotThrow()
        {
            WithLoadedMetaDraw((logic, psm, outputDirectory) =>
            {
                Assert.DoesNotThrow(() => logic.ExportSequenceCoverage(
                    new Canvas(),
                    new Canvas(),
                    outputDirectory,
                    psm));
            });
        }

        [Test]
        public void ExportAnnotatedSequence_WhenAnnotationCanvasHasNoMeasuredHeight_DoesNotThrow()
        {
            WithLoadedMetaDraw((logic, psm, outputDirectory) =>
            {
                var annotationCanvas = new Canvas();
                var legend = new Border();

                Assert.DoesNotThrow(() => logic.ExportAnnotatedSequence(
                    annotationCanvas,
                    legend,
                    psm,
                    outputDirectory,
                    200));
            });
        }

        [Test]
        public void CombineBitmap_PreservesTransparencyAndAbsoluteCoordinates()
        {
            var source = new System.Drawing.Bitmap(
                100,
                100,
                System.Drawing.Imaging.PixelFormat.Format32bppArgb);
            source.SetPixel(20, 30, System.Drawing.Color.Red);

            var combined = MetaDrawLogic.CombineBitmap(
                new List<System.Drawing.Bitmap> { source },
                new List<Point> { new(7, 11) },
                true);

            try
            {
                Assert.That(combined.GetPixel(0, 0).A, Is.EqualTo(0));
                Assert.That(combined.GetPixel(27, 41).ToArgb(), Is.EqualTo(System.Drawing.Color.Red.ToArgb()));
            }
            finally
            {
                combined.Dispose();
            }
        }

        private static void WithLoadedMetaDraw(Action<MetaDrawLogic, SpectrumMatchFromTsv, string> test)
        {
            string outputDirectory = Path.Combine(
                TestContext.CurrentContext.TestDirectory,
                "MetaDrawExportRegression_" + Guid.NewGuid().ToString("N"));
            string databasePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta");
            string spectraPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "SmallCalibratible_Yeast.mzML");

            Directory.CreateDirectory(outputDirectory);

            var searchTask = new TaskLayer.SearchTask();
            searchTask.RunTask(outputDirectory, new List<DbForTask> { new(databasePath, false) }, new List<string> { spectraPath }, "");

            string psmPath = Path.Combine(outputDirectory, "AllPSMs.psmtsv");
            var logic = new MetaDrawLogic();
            logic.SpectraFilePaths.Add(spectraPath);
            logic.SpectralMatchResultFilePaths.Add(psmPath);

            try
            {
                var errors = logic.LoadFiles(true, true);
                Assert.That(errors, Is.Empty);
                Assert.That(logic.FilteredListOfPsms, Is.Not.Empty);
                test(logic, logic.FilteredListOfPsms[0], outputDirectory);
            }
            finally
            {
                logic.CleanUpResources();
                if (Directory.Exists(outputDirectory))
                    Directory.Delete(outputDirectory, true);
            }
        }

    }
}
