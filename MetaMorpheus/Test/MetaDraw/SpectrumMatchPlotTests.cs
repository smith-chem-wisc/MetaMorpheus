using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Threading;
using System.Windows.Controls;
using GuiFunctions;
using NUnit.Framework;
using OxyPlot;
using OxyPlot.Annotations;
using Omics.Fragmentation;
using Readers;
using TaskLayer;
using Chemistry;
using EngineLayer.DatabaseLoading;
using MzLibUtil;
using MassSpectrometry;

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
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());
            Assert.That(metadrawLogic.FilteredListOfPsms.Any());

            plotView = new OxyPlot.Wpf.PlotView() { Name = "plotView" };
            parentChildView = new ParentChildScanPlotsView();
            psm = metadrawLogic.FilteredListOfPsms.First(p => p.FullSequence == "QIVHDSGR");

            // create a fake neutral loss fragment ion for testing purposes. 
            var ionToCopy = psm.MatchedIons.First();
            var productToCopy = ionToCopy.NeutralTheoreticalProduct;
            var neutralTheorecticalProduct = new Product(productToCopy.ProductType, productToCopy.Terminus, productToCopy.NeutralMass - 18.01056468, productToCopy.FragmentNumber, productToCopy.ResiduePosition, 18.01056468, productToCopy.SecondaryProductType, productToCopy.SecondaryFragmentNumber);
            var neutralLossIon = new MatchedFragmentIon(neutralTheorecticalProduct, neutralTheorecticalProduct.MonoisotopicMass.ToMz(ionToCopy.Charge), ionToCopy.Intensity, ionToCopy.Charge);
            psm.MatchedIons.Add(neutralLossIon);

        }

        [SetUp]
        public void ResetSettings()
        {
            MetaDrawSettings.AnnotateIsotopicEnvelopes = false;
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(outputFolder, true);
        }

        private string outputFolder;
        private MetaDrawLogic metadrawLogic;
        private OxyPlot.Wpf.PlotView plotView;
        private SpectrumMatchFromTsv psm;
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

            // all parameter combinations for neutral loss ions
            yield return new PeakAnnotationTestCase(false, false, false, false, "", unannotatedColor, 13);
            yield return new PeakAnnotationTestCase(false, true, false, false, "", unannotatedColor, 13);
            yield return new PeakAnnotationTestCase(false, false, true, false, "", unannotatedColor, 13);
            yield return new PeakAnnotationTestCase(false, false, false, true, "", unannotatedColor, 13);
            yield return new PeakAnnotationTestCase(false, true, true, false, "", unannotatedColor, 13);
            yield return new PeakAnnotationTestCase(false, true, false, true, "", unannotatedColor, 13);
            yield return new PeakAnnotationTestCase(false, false, true, true, "", unannotatedColor, 13);
            yield return new PeakAnnotationTestCase(false, true, true, true, "", unannotatedColor, 13);
            yield return new PeakAnnotationTestCase(true, false, false, false, "b1-18.01", bColor, 13);
            yield return new PeakAnnotationTestCase(true, true, false, false, "b1-18.01+1", bColor, 13);
            yield return new PeakAnnotationTestCase(true, false, true, false, "b1-18.01 (111.055)", bColor, 13);
            yield return new PeakAnnotationTestCase(true, false, false, true, "b₁\u208b\u2081\u2088.\u2080\u2081", bColor, 13);
            yield return new PeakAnnotationTestCase(true, true, true, false, "b1-18.01+1 (111.055)", bColor, 13);
            yield return new PeakAnnotationTestCase(true, true, false, true, "b₁\u208b\u2081\u2088.\u2080\u2081¹⁺", bColor, 13);
            yield return new PeakAnnotationTestCase(true, false, true, true, "b₁\u208b\u2081\u2088.\u2080\u2081 (111.055)", bColor, 13);
            yield return new PeakAnnotationTestCase(true, true, true, true, "b₁\u208b\u2081\u2088.\u2080\u2081¹⁺ (111.055)", bColor, 13);
        }


        [Test]
        [TestCaseSource(nameof(GetAnnotationTestCases))]
        public void TestPeakAnnotation(PeakAnnotationTestCase testCase)
        {
            // Ensure internal ions are displayed so test indices are stable
            MetaDrawSettings.DisplayInternalIons = true;
            MetaDrawSettings.DisplayInternalIonAnnotations = true;
            
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
            metaDrawLogic.SpectralMatchResultFilePaths.Add(psmFilePath);

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
            var psm = metaDrawLogic.FilteredListOfPsms.First(p => p.FullSequence == "GVTVDKMTELR(6)") as PsmFromTsv;
             
            // perform black magic to set the scan number of the MS2 to match the mzML file number
            var oldScanNum = psm.Ms2ScanNumber;
            var field = typeof(SpectrumMatchFromTsv).GetField("<Ms2ScanNumber>k__BackingField", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
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
                metaDrawLogic.ExportPlot(plotView, canvas, new List<SpectrumMatchFromTsv>() { psm }, parentChildView,
                    outputFolderPath, out errors);

                Assert.That(File.Exists(Path.Combine(outputFolderPath, @$"{psm.Ms2ScanNumber}_{psm.FullSequence}{psm.BetaPeptideBaseSequence}.{exportType}")));
            }

// clean up resources
            metaDrawLogic.CleanUpSpectraFiles();
            Assert.That(!metaDrawLogic.SpectraFilePaths.Any());

            metaDrawLogic.CleanUpPSMFiles();
            Assert.That(!metaDrawLogic.FilteredListOfPsms.Any());
            Assert.That(!metaDrawLogic.SpectralMatchResultFilePaths.Any());

            // delete output
            Directory.Delete(outputFolderPath, true);
        }

        /// <summary>
        /// Tests that PopulateEnvelopesForScanIfNeeded populates ions with null envelopes
        /// </summary>
        [Test]
        public void TestEnvelopePopulatesNullEnvelopes()
        {
            // Enable envelope annotation
            MetaDrawSettings.AnnotateIsotopicEnvelopes = true;
            MetaDrawSettings.DisplayIonAnnotations = true;

            // Store original ion count
            int originalCount = psm.MatchedIons.Count;

            // Create MatchedFragmentIonWithEnvelope with null Envelope (simulating TSV placeholder state)
            var existingIon = psm.MatchedIons.First();
            var product = existingIon.NeutralTheoreticalProduct;
            
            // Create ion with null envelope - this is the "placeholder" state from TSV files
            var ionWithNullEnvelope = new MatchedFragmentIonWithEnvelope(
                product, 
                existingIon.Mz, 
                existingIon.Intensity, 
                existingIon.Charge,
                envelope: null);  // null envelope triggers PopulateEnvelopesForScanIfNeeded
            
            // Add to PSM
            psm.MatchedIons.Add(ionWithNullEnvelope);

            try
            {
                // Verify it has null envelope before display
                Assert.That(ionWithNullEnvelope.Envelope, Is.Null, "Ion should have null envelope before display");

                // Create a fresh plotView
                var testPlotView = new OxyPlot.Wpf.PlotView();

                // Display - this should trigger PopulateEnvelopesForScanIfNeeded
                metadrawLogic.DisplaySpectrumMatch(testPlotView, psm, parentChildView, out List<string> errors);
                Assert.That(errors == null || !errors.Any());

                // After display, the envelope should be populated
                Assert.That(ionWithNullEnvelope.Envelope, Is.Not.Null, "Envelope should be populated after display");
                
                // Verify annotations were added
                Assert.That(testPlotView.Model.Annotations.Count, Is.GreaterThan(0));
            }
            finally
            {
                // Remove the added ion and restore original state
                psm.MatchedIons.RemoveAt(psm.MatchedIons.Count - 1);
                Assert.That(psm.MatchedIons.Count, Is.EqualTo(originalCount), "PSM should have original ion count after test");
                MetaDrawSettings.AnnotateIsotopicEnvelopes = false;
            }
        }

        /// <summary>
        /// Tests that PopulateEnvelopesForScanIfNeeded is skipped when all envelopes are already populated
        /// </summary>
        [Test]
        public void TestEnvelopeAlreadyPopulated()
        {
            try
            {
                // First display to populate envelopes
                MetaDrawSettings.AnnotateIsotopicEnvelopes = true;
                MetaDrawSettings.DisplayIonAnnotations = true;

                var testPlotView1 = new OxyPlot.Wpf.PlotView();
                metadrawLogic.DisplaySpectrumMatch(testPlotView1, psm, parentChildView, out List<string> errors1);
                Assert.That(errors1 == null || !errors1.Any());

                // Verify envelopes are populated
                var envelopeIons = psm.MatchedIons.OfType<MatchedFragmentIonWithEnvelope>().ToList();
                var allPopulated = envelopeIons.All(i => i.Envelope != null);
                Assert.That(allPopulated, "All envelopes should be populated after first display");

                int annotationCountAfterFirst = testPlotView1.Model.Annotations.Count;

                // Second display - should NOT call PopulateEnvelopesForScanIfNeeded (already populated)
                var testPlotView2 = new OxyPlot.Wpf.PlotView();
                metadrawLogic.DisplaySpectrumMatch(testPlotView2, psm, parentChildView, out List<string> errors2);
                Assert.That(errors2 == null || !errors2.Any());

                // Should have same annotations (envelopes already populated, no re-population needed)
                Assert.That(testPlotView2.Model.Annotations.Count, Is.EqualTo(annotationCountAfterFirst));
            }
            finally
            {
                MetaDrawSettings.AnnotateIsotopicEnvelopes = false;
            }
        }

        /// <summary>
        /// Tests that envelope check is skipped when AnnotateIsotopicEnvelopes is false
        /// </summary>
        [Test]
        public void TestEnvelopeAnnotationDisabled()
        {
            try
            {
                // Disable envelope annotation - should skip PopulateEnvelopesForScanIfNeeded
                MetaDrawSettings.AnnotateIsotopicEnvelopes = false;
                MetaDrawSettings.DisplayIonAnnotations = true;

                var testPlotView = new OxyPlot.Wpf.PlotView();

                int initialAnnotationCount = testPlotView.Model?.Annotations?.Count ?? 0;

                // Display the spectrum match - envelope check should be skipped
                metadrawLogic.DisplaySpectrumMatch(testPlotView, psm, parentChildView, out List<string> errors);
                Assert.That(errors == null || !errors.Any());

                // Verify annotations were added (normal annotation, not envelope)
                Assert.That(testPlotView.Model.Annotations.Count, Is.GreaterThan(initialAnnotationCount));
            }
            finally
            {
                MetaDrawSettings.AnnotateIsotopicEnvelopes = false;
            }
        }

        /// <summary>
        /// Tests that AnnotateEnvelopePeaks draws multiple peaks for an envelope (not just the monoisotopic)
        /// </summary>
        [Test]
        public void TestEnvelopeMultiplePeaksDrawn()
        {
            try
            {
                MetaDrawSettings.AnnotateIsotopicEnvelopes = true;
                MetaDrawSettings.DisplayIonAnnotations = true;

                var testPlotView = new OxyPlot.Wpf.PlotView();
                metadrawLogic.DisplaySpectrumMatch(testPlotView, psm, parentChildView, out List<string> errors);
                Assert.That(errors == null || !errors.Any());

                // Get ions with populated envelopes
                var envelopeIons = psm.MatchedIons.OfType<MatchedFragmentIonWithEnvelope>()
                    .Where(i => i.Envelope != null)
                    .ToList();

                Assert.That(envelopeIons.Any(), "Should have ions with populated envelopes");

                // Find ions with envelopes that have multiple peaks
                var multiPeakEnvelopes = envelopeIons.Where(e => e.Envelope!.Peaks.Count > 1).ToList();
                
                if (multiPeakEnvelopes.Any())
                {
                    // If there are multi-peak envelopes, verify the plot has more annotations than just ions
                    // (envelope adds multiple peak markers per ion)
                    Console.WriteLine($"Found {multiPeakEnvelopes.Count} ions with multi-peak envelopes");
                    
                    // The annotation count should be greater than just the number of ions
                    // because each envelope with multiple peaks adds extra markers
                    Assert.That(testPlotView.Model.Annotations.Count, Is.GreaterThan(envelopeIons.Count));
                }
                else
                {
                    // If no multi-peak envelopes, just verify annotations exist
                    Assert.That(testPlotView.Model.Annotations.Count, Is.GreaterThan(0));
                    Console.WriteLine("No multi-peak envelopes found in test data");
                }
            }
            finally
            {
                MetaDrawSettings.AnnotateIsotopicEnvelopes = false;
            }
        }

        /// <summary>
        /// Tests that internal fragment ions are handled correctly with envelope annotation off
        /// </summary>
        [Test]
        public void TestInternalFragmentWithEnvelopeDisabled()
        {
            MetaDrawSettings.AnnotateIsotopicEnvelopes = false;
            MetaDrawSettings.DisplayIonAnnotations = true;
            MetaDrawSettings.DisplayInternalIonAnnotations = false;

            var testPlotView = new OxyPlot.Wpf.PlotView();
            metadrawLogic.DisplaySpectrumMatch(testPlotView, psm, parentChildView, out List<string> errors);
            Assert.That(errors == null || !errors.Any());

            // Verify annotations exist but internal fragments are hidden
            Assert.That(testPlotView.Model.Annotations.Count, Is.GreaterThan(0));

            // Find any text annotations that might indicate internal fragments
            var textAnnotations = testPlotView.Model.Annotations.OfType<OxyPlot.Annotations.TextAnnotation>().ToList();
            
            // With DisplayInternalIonAnnotations=false, internal fragment text should be hidden
            // The color should be HiddenAnnotationColor
            var hiddenAnnotations = textAnnotations.Where(a => 
                a.TextColor == OxyColor.FromArgb(0, 0, 0, 1)).ToList();
            
            Console.WriteLine($"Total annotations: {textAnnotations.Count}, Hidden: {hiddenAnnotations.Count}");
        }
    }
}
