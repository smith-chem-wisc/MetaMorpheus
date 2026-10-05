using EngineLayer;
using EngineLayer.DatabaseLoading;
using GuiFunctions;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
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

        [TestCase("Png")]
        [TestCase("Tiff")]
        [TestCase("Wmf")]
        [TestCase("Pdf")]
        public void CombineBitmap_AlphaCapableFormats_KeepTransparentBackground(string exportType)
        {
            WithExportType(exportType, () =>
            {
                var source = new System.Drawing.Bitmap(
                    50,
                    50,
                    System.Drawing.Imaging.PixelFormat.Format32bppArgb);
                source.SetPixel(10, 10, System.Drawing.Color.Red);

                var combined = MetaDrawLogic.CombineBitmap(
                    new List<System.Drawing.Bitmap> { source },
                    new List<Point> { new(0, 0) },
                    true);

                try
                {
                    Assert.That(combined.GetPixel(0, 0).A, Is.EqualTo(0),
                        $"Uncovered background should stay transparent for {exportType}.");
                    Assert.That(combined.GetPixel(10, 10).ToArgb(), Is.EqualTo(System.Drawing.Color.Red.ToArgb()),
                        $"Drawn pixels should still appear for {exportType}.");
                }
                finally
                {
                    combined.Dispose();
                }
            });
        }

        [TestCase("Jpeg")]
        [TestCase("Bmp")]
        public void CombineBitmap_NonAlphaFormats_PaintOpaqueWhiteBackground(string exportType)
        {
            WithExportType(exportType, () =>
            {
                var source = new System.Drawing.Bitmap(
                    50,
                    50,
                    System.Drawing.Imaging.PixelFormat.Format32bppArgb);
                source.SetPixel(10, 10, System.Drawing.Color.Red);

                var combined = MetaDrawLogic.CombineBitmap(
                    new List<System.Drawing.Bitmap> { source },
                    new List<Point> { new(0, 0) },
                    true);

                try
                {
                    var uncovered = combined.GetPixel(0, 0);
                    Assert.That(uncovered.A, Is.EqualTo(255),
                        $"Uncovered background alpha should be opaque for {exportType}.");
                    Assert.That(uncovered.ToArgb(), Is.EqualTo(System.Drawing.Color.White.ToArgb()),
                        $"Uncovered background should be white for {exportType}.");
                    Assert.That(combined.GetPixel(10, 10).ToArgb(), Is.EqualTo(System.Drawing.Color.Red.ToArgb()),
                        $"Drawn pixels should still appear for {exportType}.");
                }
                finally
                {
                    combined.Dispose();
                }
            });
        }

        [Test]
        public void CombineBitmap_ExplicitBackground_OverridesDefault()
        {
            var source = new System.Drawing.Bitmap(
                20,
                20,
                System.Drawing.Imaging.PixelFormat.Format32bppArgb);

            var combined = MetaDrawLogic.CombineBitmap(
                new List<System.Drawing.Bitmap> { source },
                new List<Point> { new(0, 0) },
                true,
                System.Drawing.Color.Lime);

            try
            {
                Assert.That(combined.GetPixel(0, 0).ToArgb(), Is.EqualTo(System.Drawing.Color.Lime.ToArgb()));
            }
            finally
            {
                combined.Dispose();
            }
        }

        [TestCase("Jpeg")]
        [TestCase("Bmp")]
        public void ExportSequenceCoverage_NonAlphaFormats_ProduceOpaqueWhiteFile(string exportType)
        {
            WithExportType(exportType, () =>
            {
                WithLoadedMetaDraw((logic, psm, outputDirectory) =>
                {
                    var textCanvas = new Canvas();
                    var mapCanvas = new Canvas();
                    logic.DisplaySequences(mapCanvas, null, null, psm);

                    Assert.DoesNotThrow(() => logic.ExportSequenceCoverage(
                        textCanvas,
                        mapCanvas,
                        outputDirectory,
                        psm));

                    string filePath = Directory.GetFiles(outputDirectory, "*_SequenceCoverage.*").Single();
                    using var bitmap = new System.Drawing.Bitmap(filePath);
                    Assert.That(bitmap.GetPixel(0, 0).A, Is.EqualTo(255),
                        $"{exportType} output should not contain transparent pixels.");
                });
            });
        }

        private static void WithExportType(string exportType, Action body)
        {
            string originalExportType = MetaDrawSettings.ExportType;
            try
            {
                MetaDrawSettings.ExportType = exportType;
                body();
            }
            finally
            {
                MetaDrawSettings.ExportType = originalExportType;
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
