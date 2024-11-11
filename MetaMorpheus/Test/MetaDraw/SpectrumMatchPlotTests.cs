using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Windows.Controls;
using EngineLayer;
using GuiFunctions;
using NUnit.Framework;
using OxyPlot;
using OxyPlot.Annotations;
using Omics.Fragmentation;
using Readers;
using TaskLayer;

namespace Test.MetaDraw
{
    /// <summary>
    /// Class for testing different settings and how they affect the displayed elements
    /// The set up will instantiate everything needed to plot a simple psm, then test cases can
    /// be written to test different settings
    /// </summary>
    [TestFixture, Apartment(ApartmentState.STA)]
    public class SpectrumMatchPlotTests
    {
        /// <summary>
        /// Run search task and prepare everything needed to plot a simple psm
        /// </summary>
        [OneTimeSetUp]
        public void Setup()
        {
            
            outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"MetaDraw_PeakAnnotaitonTest");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");

            Directory.CreateDirectory(outputFolder);
             
            // run search task
            var searchtask = new SearchTask()
            {
                SearchParameters = new SearchParameters()
                {
                    MinAllowedInternalFragmentLength = 4,
                }
            };
            searchtask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, new List<string> { spectraFile }, "");

            var psmFile = Path.Combine(outputFolder, @"AllPSMs.psmtsv");

            // load results into metadraw
            metadrawLogic = new MetaDrawLogic();
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());
            Assert.That(metadrawLogic.FilteredListOfPsms.Any());

            plotView = new OxyPlot.Wpf.PlotView() { Name = "plotView" };
            parentChildView = new ParentChildScanPlotsView();
            psm = metadrawLogic.FilteredListOfPsms.First(p => p.FullSequence == "QIVHDSGR");
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(outputFolder, true);
        }

        private string outputFolder;
        private MetaDrawLogic metadrawLogic;
        private OxyPlot.Wpf.PlotView plotView;
        private PsmFromTsv psm;
        private ParentChildScanPlotsView parentChildView;
        public record PeakAnnotationTestCase(bool AnnotatePeaks, bool AnnotateCharges,
            bool AnnotateMass, bool SubAndSuper, string ExpectedAnnotation, OxyColor ExpectedColor, int FragmentIndex);

        public static IEnumerable<PeakAnnotationTestCase> GetAnnotationTestCases()
        {
            OxyColor unannotatedColor = OxyColor.FromArgb(0, 0, 0, 1);
            OxyColor bColor = OxyColors.Blue;
            OxyColor yColor = OxyColors.Red;
            OxyColor internalColor = OxyColors.Purple;

            // all parameter combinations for b ions
            yield return new PeakAnnotationTestCase(false, false, false, false, "", unannotatedColor, 0);
            yield return new PeakAnnotationTestCase(false, true, false, false, "", unannotatedColor, 0);
            yield return new PeakAnnotationTestCase(false, false, true, false, "", unannotatedColor, 0);
            yield return new PeakAnnotationTestCase(false, false, false, true, "", unannotatedColor, 0);
            yield return new PeakAnnotationTestCase(false, true, true, false, "", unannotatedColor, 0);
            yield return new PeakAnnotationTestCase(false, true, false, true, "", unannotatedColor, 0);
            yield return new PeakAnnotationTestCase(false, false, true, true, "", unannotatedColor, 0);
            yield return new PeakAnnotationTestCase(false, true, true, true, "", unannotatedColor, 0);
            yield return new PeakAnnotationTestCase(true, false, false, false, "b1", bColor, 0);
            yield return new PeakAnnotationTestCase(true, true, false, false, "b1+1", bColor, 0);
            yield return new PeakAnnotationTestCase(true, false, true, false, "b1 (129.066)", bColor, 0);
            yield return new PeakAnnotationTestCase(true, false, false, true, "b₁", bColor, 0);
            yield return new PeakAnnotationTestCase(true, true, true, false, "b1+1 (129.066)", bColor, 0);
            yield return new PeakAnnotationTestCase(true, true, false, true, "b₁¹⁺", bColor, 0);
            yield return new PeakAnnotationTestCase(true, false, true, true, "b₁ (129.066)", bColor, 0);
            yield return new PeakAnnotationTestCase(true, true, true, true, "b₁¹⁺ (129.066)", bColor, 0);

            // all parameter combinations for y ions
            yield return new PeakAnnotationTestCase(false, false, false, false, "", unannotatedColor, 9);
            yield return new PeakAnnotationTestCase(false, true, false, false, "", unannotatedColor, 9);
            yield return new PeakAnnotationTestCase(false, false, true, false, "", unannotatedColor, 9);
            yield return new PeakAnnotationTestCase(false, false, false, true, "", unannotatedColor, 9);
            yield return new PeakAnnotationTestCase(false, true, true, false, "", unannotatedColor, 9);
            yield return new PeakAnnotationTestCase(false, true, false, true, "", unannotatedColor, 9);
            yield return new PeakAnnotationTestCase(false, false, true, true, "", unannotatedColor, 9);
            yield return new PeakAnnotationTestCase(false, true, true, true, "", unannotatedColor, 9);
            yield return new PeakAnnotationTestCase(true, false, false, false, "y7", yColor, 9);
            yield return new PeakAnnotationTestCase(true, true, false, false, "y7+1", yColor, 9);
            yield return new PeakAnnotationTestCase(true, false, true, false, "y7 (783.411)", yColor, 9);
            yield return new PeakAnnotationTestCase(true, false, false, true, "y₇", yColor, 9);
            yield return new PeakAnnotationTestCase(true, true, true, false, "y7+1 (783.411)", yColor, 9);
            yield return new PeakAnnotationTestCase(true, true, false, true, "y₇¹⁺", yColor, 9);
            yield return new PeakAnnotationTestCase(true, false, true, true, "y₇ (783.411)", yColor, 9);
            yield return new PeakAnnotationTestCase(true, true, true, true, "y₇¹⁺ (783.411)", yColor, 9);

            // all parameter combinations for internal ions
            yield return new PeakAnnotationTestCase(false, false, false, false, "", unannotatedColor, 12);
            yield return new PeakAnnotationTestCase(false, true, false, false, "", unannotatedColor, 12);
            yield return new PeakAnnotationTestCase(false, false, true, false, "", unannotatedColor, 12);
            yield return new PeakAnnotationTestCase(false, false, false, true, "", unannotatedColor, 12);
            yield return new PeakAnnotationTestCase(false, true, true, false, "", unannotatedColor, 12);
            yield return new PeakAnnotationTestCase(false, true, false, true, "", unannotatedColor, 12);
            yield return new PeakAnnotationTestCase(false, false, true, true, "", unannotatedColor, 12);
            yield return new PeakAnnotationTestCase(false, true, true, true, "", unannotatedColor, 12);
            yield return new PeakAnnotationTestCase(true, false, false, false, "yIb[4-7]", internalColor, 12);
            yield return new PeakAnnotationTestCase(true, true, false, false, "yIb[4-7]+1", internalColor, 12);
            yield return new PeakAnnotationTestCase(true, false, true, false, "yIb[4-7] (397.145)", internalColor, 12);
            yield return new PeakAnnotationTestCase(true, false, false, true, "yIb₄₋₇", internalColor, 12);
            yield return new PeakAnnotationTestCase(true, true, true, false, "yIb[4-7]+1 (397.145)", internalColor, 12);
            yield return new PeakAnnotationTestCase(true, true, false, true, "yIb₄₋₇¹⁺", internalColor, 12);
            yield return new PeakAnnotationTestCase(true, false, true, true, "yIb₄₋₇ (397.145)", internalColor, 12);
            yield return new PeakAnnotationTestCase(true, true, true, true, "yIb₄₋₇¹⁺ (397.145)", internalColor, 12);
        }


        [Test]
        [TestCaseSource(nameof(GetAnnotationTestCases))]
        public void TestPeakAnnotation(PeakAnnotationTestCase testCase)
        {
            // set parameters for test case
            MetaDrawSettings.DisplayIonAnnotations = testCase.AnnotatePeaks;
            MetaDrawSettings.AnnotateCharges = testCase.AnnotateCharges;
            MetaDrawSettings.SubAndSuperScriptIons = testCase.SubAndSuper;
            MetaDrawSettings.AnnotateMzValues = testCase.AnnotateMass;

            metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out List<string> errors);
            Assert.That(errors == null || !errors.Any());

            var annotation = ((TextAnnotation)plotView.Model.Annotations[testCase.FragmentIndex]).Text;
            Assert.That(annotation, Is.EqualTo(testCase.ExpectedAnnotation));
        }

        [Test]
        [TestCaseSource(nameof(GetAnnotationTestCases))]
        public void TestPeakColor(PeakAnnotationTestCase testCase)
        {
            MetaDrawSettings.DisplayInternalIons = true;
            MetaDrawSettings.DisplayInternalIonAnnotations = true;
            MetaDrawSettings.DisplayIonAnnotations = testCase.AnnotatePeaks;
            MetaDrawSettings.AnnotateCharges = testCase.AnnotateCharges;
            MetaDrawSettings.SubAndSuperScriptIons = testCase.SubAndSuper;
            MetaDrawSettings.AnnotateMzValues = testCase.AnnotateMass;

            metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out List<string> errors);
            Assert.That(errors == null || !errors.Any());

            var annotation = ((TextAnnotation)plotView.Model.Annotations[testCase.FragmentIndex]);
            Assert.That(annotation.TextColor, Is.EqualTo(testCase.ExpectedColor));
        }

        [Test]
        public static void TestCrossLinkSpectrumMatchPlot()
        { 
            // set up file paths
            var outputFolderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestCrossLinkSpectrumMatchPlot");
            var psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "XlTestData", "XL_Interlinks.tsv");
            var dataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "XlTestData", "2017-11-21_XL_DSSO_Ribosome_RT60min_28800-28898.mzML");

            Directory.CreateDirectory(outputFolderPath);

            // load in files
            MetaDrawLogic metaDrawLogic = new MetaDrawLogic();
            metaDrawLogic.SpectraFilePaths.Add(dataFilePath);
            metaDrawLogic.PsmResultFilePaths.Add(psmFilePath);

            var errors = metaDrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());
            Assert.That(metaDrawLogic.FilteredListOfPsms.Any());

            // load in gui components
            var plotView = new OxyPlot.Wpf.PlotView() { Name = "plotView" };
            var canvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var stationaryCanvas = new Canvas();
            var sequenceAnnotationCanvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            var psm = metaDrawLogic.FilteredListOfPsms.First(p => p.FullSequence == "GVTVDKMTELR(6)");
             
            // perform black magic to set the scan number of the MS2 to match the mzML file number
            var oldScanNum = psm.Ms2ScanNumber;
            var field = typeof(PsmFromTsv).GetField("<Ms2ScanNumber>k__BackingField", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            field?.SetValue(psm, 28819);

            // display psm and check display has correct number of annotations
            metaDrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, psm);
            int alphaPeptideAnnotations = psm.BaseSeq.Length + psm.MatchedIons.Count(p => p.NeutralTheoreticalProduct.ProductType is ProductType.b or ProductType.y);
            int betaPeptideAnnotations = psm.BetaPeptideBaseSequence.Length + psm.BetaPeptideMatchedIons.Count(p => p.NeutralTheoreticalProduct.ProductType is ProductType.b or ProductType.y);
            int crossLinkerAnnotations = 1;
            int sequenceDisplayAnnotationCount =
                alphaPeptideAnnotations + betaPeptideAnnotations + crossLinkerAnnotations;
            Assert.That(stationaryCanvas.Children.Count, Is.EqualTo(sequenceDisplayAnnotationCount));
            Assert.That(scrollableCanvas.Children.Count, Is.EqualTo(alphaPeptideAnnotations));

            var alphaPeptideSpectralMatchAnnotationCount = psm.MatchedIons.Count;
            var betaPeptideSpectralMatchAnnotationCount = psm.BetaPeptideMatchedIons.Count;
            var csmPlotIonAnnotationCount = alphaPeptideSpectralMatchAnnotationCount + betaPeptideSpectralMatchAnnotationCount;
            metaDrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out errors);
            Assert.That(plotView.Model.Annotations.Count, Is.EqualTo(csmPlotIonAnnotationCount + 1)); // Plus One for the annotation text describing the csm

            var scan = MsDataFileReader.GetDataFile(dataFilePath).LoadAllStaticData().GetOneBasedScan(oldScanNum);
            var csmPlot = new CrosslinkSpectrumMatchPlot(plotView, psm, scan, stationaryCanvas, false);
            Assert.That(csmPlot.Model.Annotations.Count, Is.EqualTo(csmPlotIonAnnotationCount));

            // test each export type
            foreach (var exportType in MetaDrawSettings.ExportTypes)
            {
                MetaDrawSettings.ExportType = exportType;
                metaDrawLogic.ExportPlot(plotView, canvas, new List<PsmFromTsv>() { psm }, parentChildView,
                    outputFolderPath, out errors);

                Assert.That(File.Exists(Path.Combine(outputFolderPath, @$"{psm.Ms2ScanNumber}_{psm.FullSequence}{psm.BetaPeptideBaseSequence}.{exportType}")));
            }

            // clean up resources
            metaDrawLogic.CleanUpSpectraFiles();
            Assert.That(!metaDrawLogic.SpectraFilePaths.Any());

            metaDrawLogic.CleanUpPSMFiles();
            Assert.That(!metaDrawLogic.FilteredListOfPsms.Any());
            Assert.That(!metaDrawLogic.PsmResultFilePaths.Any());

            // delete output
            Directory.Delete(outputFolderPath, true);
        }
    }
}
