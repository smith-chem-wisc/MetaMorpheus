using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Controls;
using EngineLayer;
using GuiFunctions;
using NUnit.Framework;
using OxyPlot.Annotations;
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

        [OneTimeSetUp]
        public void Setup()
        {
            outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"MetaDraw_PeakAnnotaitonTest");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");

            Directory.CreateDirectory(outputFolder);

            // run search task
            var searchtask = new SearchTask();
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

        private string outputFolder;
        private MetaDrawLogic metadrawLogic;
        private OxyPlot.Wpf.PlotView plotView;
        private PsmFromTsv psm;
        private ParentChildScanPlotsView parentChildView;

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(outputFolder, true);
        }

        [Test]
        [TestCase(false, false, false, false, "")]
        [TestCase(false, true, false, false, "")]
        [TestCase(false, false, true, false, "")]
        [TestCase(false, false, false, true, "")]
        [TestCase(false, true, true, false, "")]
        [TestCase(false, true, false, true, "")]
        [TestCase(false, false, true, true, "")]
        [TestCase(false, true, true, true, "")]
        [TestCase(true, false, false, false, "b1")]
        [TestCase(true, true, false, false, "b1+1")]
        [TestCase(true, false, true, false, "b1 (129.066)")]
        [TestCase(true, false, false, true, "b₁")]
        [TestCase(true, true, true, false, "b1+1 (129.066)")]
        [TestCase(true, true, false, true, "b₁¹⁺")]
        [TestCase(true, false, true, true, "b₁ (129.066)")]
        [TestCase(true, true, true, true, "b₁¹⁺ (129.066)")]
        public void TestPeakAnnotation(bool annotatePeaks, bool annotateCharges, bool annotateMass, bool subAndSuper,
            string expected)

        {
            // set parameters for test case
            MetaDrawSettings.DisplayIonAnnotations = annotatePeaks;
            MetaDrawSettings.AnnotateCharges = annotateCharges;
            MetaDrawSettings.SubAndSuperScriptIons = subAndSuper;
            MetaDrawSettings.AnnotateMzValues = annotateMass;
            metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out List<string> errors);
            Assert.That(errors == null || !errors.Any());

            var annotation = ((TextAnnotation)plotView.Model.Annotations.First()).Text;
            Assert.That(annotation, Is.EqualTo(expected));

        }
    }
}
