using Chemistry;
using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using OxyPlot.Series;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Windows.Controls;
using TaskLayer;

namespace Test
{
    [TestFixture, Apartment(ApartmentState.STA)]
    [TestFixture]
    public class MetaDrawTest
    {
        [Test]
        public static void TestMetaDrawReadPsmFile()
        {
            SearchTask searchTask = new SearchTask();

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMetaDrawReadPsmFile");

            DbForTask db = new DbForTask(myDatabase, false);
            Directory.CreateDirectory(folderPath);

            searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "metadraw");
            string psmFile = Directory.GetFiles(folderPath).First(f => f.Contains("AllPSMs.psmtsv"));

            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);

            Assert.AreEqual(10, parsedPsms.Count);
            Assert.AreEqual(0, warnings.Count);

            Directory.Delete(folderPath, true);
        }

        [Test]
        public static void TestParenthesesRemovalForSilac()
        {
            string baseSequence = "ASDF(+8.01)ASDF";
            string cleanedSequence = PsmFromTsv.RemoveParentheses(baseSequence);
            Assert.IsTrue(cleanedSequence.Equals("ASDFASDF"));
        }

        [Test]
        public static void TestMetaDrawReadCrossPsmFile()
        {
            XLSearchTask searchTask = new XLSearchTask();
            searchTask.XlSearchParameters.Crosslinker = GlobalVariables.Crosslinkers.ToList()[1];

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA_DSS_23747.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMetaDrawReadPsmFile");

            DbForTask db = new DbForTask(myDatabase, false);
            Directory.CreateDirectory(folderPath);

            searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "metadraw");

            string psmFile = Directory.GetFiles(folderPath).First(f => f.Contains("XL_Intralinks.tsv"));

            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);

            Directory.Delete(folderPath, true);
        }

        [Test]
        public static void TestFindVariantCrossingIons()
        {
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\VariantCrossTest.psmtsv");
            List<string> warnings = new List<string>();
            List<PsmFromTsv> psms;

            psms = PsmTsvReader.ReadTsv(myFile, out warnings);  // test will fail if the order of psms is changed to something other than top to bottom row

            // check that variant psm properties are being parsed correctly
            Assert.AreEqual("", psms[0].IdentifiedSequenceVariations);
            Assert.AreEqual("A147T", psms[1].IdentifiedSequenceVariations);

            Assert.AreEqual("541-541", psms[0].SpliceSites);
            Assert.AreEqual("", psms[1].SpliceSites);

            // check that the correct ions are being added to VariantCrossingIons
            List<List<string>> expected = new List<List<string>>();
            expected.Add(new List<string>() { });                           // no variant (7527)
            expected.Add(new List<string>() { "b3" });                      // b fragments before and after variant (4211)
            expected.Add(new List<string>() { "y3" });                      // y fragments before and after variant (7759)
            expected.Add(new List<string>() { "b4", "b16" });               // variant at 1st position (9221)
            expected.Add(new List<string>() { "y1", "y9" });                // variant at last position (6778)
            expected.Add(new List<string>() { "b4", "y35" });               // peptide comes from multiple proteins (8613)
            expected.Add(new List<string>() { "b1", "b10", "y1", "y6" });   // variation spans the whole peptide (8765)
            expected.Add(new List<string>() { });                           // variant before peptide (8169)
            expected.Add(new List<string>() { });                           // variant after peptide (6532)
            expected.Add(new List<string>() { "a3" });                      // a fragments before and after variant (4212)
            expected.Add(new List<string>() { "c3" });                      // c fragments before and after variant (4213)
            expected.Add(new List<string>() { "x3" });                      // x fragments before and after variant (7760)
            expected.Add(new List<string>() { "zDot3" });                   // z fragments before and after variant (7761)
            expected.Add(new List<string>() { });   // M fragment with length almost the whole peptide and variant in the middle (7762)
            expected.Add(new List<string>() { });   // D fragment with length almost the whole peptide and variant in the middle (7763)

            for (int i = 0; i < psms.Count; i++)
            {
                IEnumerable<string> actualIons = psms[i].VariantCrossingIons.Select(p => p.NeutralTheoreticalProduct.Annotation);
                foreach (string expectedIon in expected[i])
                    Assert.IsTrue(actualIons.Contains(expectedIon),
                       "VariantCrossingIons should contain ion " + expectedIon + " in file " + psms[i].FileNameWithoutExtension + ".");
                foreach (string actualIon in actualIons)
                    Assert.IsTrue(expected[i].Contains(actualIon),
                        "VariantCrossingIons should not contain ion " + actualIon + " in file " + psms[i].FileNameWithoutExtension + ".");
            }
        }

        [Test]
        public static void MetaDraw_SearchTaskTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"MetaDraw_SearchTaskTest");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");

            Directory.CreateDirectory(outputFolder);

            // run search task
            var searchtask = new SearchTask();
            searchtask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, new List<string> { spectraFile }, "");

            var psmFile = Path.Combine(outputFolder, @"AllPSMs.psmtsv");

            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());
            Assert.That(metadrawLogic.FilteredListOfPsms.Any());

            // test results filter
            MetaDrawSettings.QValueFilter = 0.01;
            MetaDrawSettings.ShowDecoys = false;
            metadrawLogic.FilterPsms();
            Assert.That(metadrawLogic.FilteredListOfPsms.All(p => p.DecoyContamTarget == "T"));
            Assert.That(metadrawLogic.FilteredListOfPsms.All(p => p.QValue <= 0.01));

            MetaDrawSettings.QValueFilter = 1.0;
            MetaDrawSettings.ShowDecoys = true;
            metadrawLogic.FilterPsms();
            Assert.That(metadrawLogic.FilteredListOfPsms.Any(p => p.DecoyContamTarget == "D"));
            Assert.That(metadrawLogic.FilteredListOfPsms.Any(p => p.QValue > 0.01));

            // test text search filter (filter by full sequence)
            string filterString = @"QIVHDSGR";
            metadrawLogic.FilterPsmsByString(filterString);

            int c = 0;
            foreach (var filteredPsm in metadrawLogic.PeptideSpectralMatchesView)
            {
                var psmObj = (PsmFromTsv)filteredPsm;
                Assert.That(psmObj.FullSequence.Contains(filterString));
                c++;
            }
            Assert.That(c > 0);

            // test text search filter (filter by MS2 scan number)
            filterString = @"120";
            metadrawLogic.FilterPsmsByString(filterString);

            c = 0;
            foreach (var filteredPsm in metadrawLogic.PeptideSpectralMatchesView)
            {
                var psmObj = (PsmFromTsv)filteredPsm;
                Assert.That(psmObj.Ms2ScanNumber.ToString().Contains(filterString));
                c++;
            }
            Assert.That(c > 0);

            // draw PSM
            var plotView = new OxyPlot.Wpf.PlotView();
            var canvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            var psm = metadrawLogic.FilteredListOfPsms.First();

            metadrawLogic.DisplaySpectrumMatch(plotView, canvas, psm, parentChildView, out errors);
            Assert.That(errors == null || !errors.Any());

            // test that plot was drawn
            var plotSeries = plotView.Model.Series;
            var series = plotSeries[0]; // the first m/z peak
            var peakPoints = ((LineSeries)series).Points;
            Assert.That(Math.Round(peakPoints[0].X, 2) == 101.07); // m/z
            Assert.That(Math.Round(peakPoints[1].X, 2) == 101.07);
            Assert.That((int)peakPoints[0].Y == 0); // intensity
            Assert.That((int)peakPoints[1].Y == 35045);

            var plotAxes = plotView.Model.Axes;
            Assert.That(plotAxes.Count == 2);

            // test that base sequence annotation was drawn
            int numAnnotatedResidues = psm.BaseSeq.Length;
            int numAnnotatedIons = psm.MatchedIons.Count;
            int numAnnotatedMods = psm.FullSequence.Count(p => p == '[');
            Assert.That(canvas.Children.Count == numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods);

            // write pdf
            var psmsToExport = metadrawLogic.FilteredListOfPsms.Where(p => p.FullSequence == "QIVHDSGR").Take(3).ToList();
            metadrawLogic.ExportToPdf(plotView, canvas, psmsToExport, parentChildView, outputFolder, out errors);

            // test that pdf exists
            Assert.That(File.Exists(Path.Combine(outputFolder, @"116_QIVHDSGR.pdf")));
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR.pdf")));
            Assert.That(File.Exists(Path.Combine(outputFolder, @"127_QIVHDSGR.pdf")));

            // test displaying a PSM with a mod
            var modPsm = metadrawLogic.FilteredListOfPsms.First(p => p.FullSequence.Contains("["));
            metadrawLogic.DisplaySpectrumMatch(plotView, canvas, modPsm, parentChildView, out errors);
            Assert.That(errors == null || !errors.Any());
            Assert.That(canvas.Children.Count == modPsm.BaseSeq.Length + modPsm.MatchedIons.Count + modPsm.FullSequence.Count(p => p == '['));

            // clean up resources
            metadrawLogic.CleanUpResources();
            Assert.That(!metadrawLogic.FilteredListOfPsms.Any());
            Assert.That(!metadrawLogic.PsmResultFilePaths.Any());
            Assert.That(!metadrawLogic.SpectraFilePaths.Any());

            // delete output
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void MetaDraw_XlSearchTaskWithChildScansTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"MetaDraw_XlSearchTaskTest");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA.fasta");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\ms2mixed_bsa_xlink.mzML");

            // run task
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.CID, ms2childScanDissociationType: DissociationType.ETD,
                trimMsMsPeaks: false);

            Directory.CreateDirectory(outputFolder);
            var xlSearchTask = new XLSearchTask() { CommonParameters = commonParameters };
            xlSearchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, new List<string> { spectraFile }, "");

            //TODO: test other files (XL_Interlinks.tsv, Deadends.tsv, Looplinks.tsv, SinglePeptides.tsv)
            var csmFile = Path.Combine(outputFolder, @"XL_Intralinks.tsv");

            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            metadrawLogic.PsmResultFilePaths.Add(csmFile);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());
            Assert.That(metadrawLogic.FilteredListOfPsms.Any());

            // test results filter
            MetaDrawSettings.QValueFilter = 0.01;
            MetaDrawSettings.ShowDecoys = false;
            metadrawLogic.FilterPsms();
            Assert.That(metadrawLogic.FilteredListOfPsms.All(p => p.DecoyContamTarget == "T"));
            Assert.That(metadrawLogic.FilteredListOfPsms.All(p => p.QValue <= 0.01));

            // test text search filter (filter by full sequence)
            string filterString = @"SLGKVGTR";
            metadrawLogic.FilterPsmsByString(filterString);

            int c = 0;
            foreach (var filteredPsm in metadrawLogic.PeptideSpectralMatchesView)
            {
                var psmObj = (PsmFromTsv)filteredPsm;
                Assert.That(psmObj.FullSequence.Contains(filterString));
                c++;
            }
            Assert.That(c > 0);

            // test text search filter (filter by MS2 scan number)
            filterString = @"2";
            metadrawLogic.FilterPsmsByString(filterString);

            c = 0;
            foreach (var filteredPsm in metadrawLogic.PeptideSpectralMatchesView)
            {
                var psmObj = (PsmFromTsv)filteredPsm;
                Assert.That(psmObj.Ms2ScanNumber.ToString().Contains(filterString));
                c++;
            }
            Assert.That(c > 0);

            // draw PSM
            var plotView = new OxyPlot.Wpf.PlotView();
            var canvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            var csm = metadrawLogic.FilteredListOfPsms.First();

            metadrawLogic.DisplaySpectrumMatch(plotView, canvas, csm, parentChildView, out errors);
            Assert.That(errors == null || !errors.Any());

            // test that plot was drawn
            var peak = (LineSeries)plotView.Model.Series[0]; // the first m/z peak
            var peakPoints = peak.Points;
            Assert.That(Math.Round(peakPoints[0].X, 2) == 142.12); // m/z
            Assert.That(Math.Round(peakPoints[1].X, 2) == 142.12);
            Assert.That((int)peakPoints[0].Y == 0); // intensity
            Assert.That((int)peakPoints[1].Y == 1114);

            var plotAxes = plotView.Model.Axes;
            Assert.That(plotAxes.Count == 2);

            // test that base sequence annotation was drawn
            Assert.That(canvas.Children.Count > 0);

            // test that the plots were drawn in the parent/child view
            Assert.That(parentChildView.Plots.Count == 2);

            // test parent scan
            var parentPlot = parentChildView.Plots[0];
            Assert.That(parentPlot.SpectrumLabel == "Scan: 2 Dissociation Type: CID MsOrder: 2 Selected Mz: 492.02 Retention Time: 23.9");
            int numAnnotatedResidues = csm.BaseSeq.Length;
            int numAnnotatedIons = csm.MatchedIons.Count(p => p.NeutralTheoreticalProduct.ProductType != ProductType.M);
            int numAnnotatedMods = csm.FullSequence.Count(p => p == '[');
            Assert.That(parentPlot.TheCanvas.Children.Count == numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods);

            peak = (LineSeries)parentPlot.Plot.Model.Series[0]; // the first m/z peak
            peakPoints = peak.Points;
            Assert.That(Math.Round(peakPoints[0].X, 2) == 142.12); // m/z
            Assert.That(Math.Round(peakPoints[1].X, 2) == 142.12);
            Assert.That((int)peakPoints[0].Y == 0); // intensity
            Assert.That((int)peakPoints[1].Y == 1114);

            // test child scan
            var childPlot = parentChildView.Plots[1];
            Assert.That(childPlot.SpectrumLabel == "Scan: 3 Dissociation Type: ETD MsOrder: 2 Selected Mz: 492.02 RetentionTime: 23.9");
            Assert.That(childPlot.TheCanvas.Children.Count > 0);
            numAnnotatedResidues = csm.BaseSeq.Length;
            numAnnotatedIons = csm.ChildScanMatchedIons[3].Concat(csm.BetaPeptideChildScanMatchedIons[3])
                .Count(p => p.NeutralTheoreticalProduct.ProductType != ProductType.M);
            numAnnotatedMods = csm.FullSequence.Count(p => p == '[');
            Assert.That(childPlot.TheCanvas.Children.Count == numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods);

            peak = (LineSeries)childPlot.Plot.Model.Series[0]; // the first m/z peak
            peakPoints = peak.Points;
            Assert.That(Math.Round(peakPoints[0].X, 2) == 122.92); // m/z
            Assert.That(Math.Round(peakPoints[1].X, 2) == 122.92);
            Assert.That((int)peakPoints[0].Y == 0); // intensity
            Assert.That((int)peakPoints[1].Y == 857);

            // write pdf
            var psmsToExport = metadrawLogic.FilteredListOfPsms.Where(p => p.FullSequence == "SLGKVGTR(4)").ToList();
            metadrawLogic.ExportToPdf(plotView, canvas, psmsToExport, parentChildView, outputFolder, out errors);

            // test that pdf exists
            Assert.That(File.Exists(Path.Combine(outputFolder, @"2_SLGKVGTR(4).pdf"))); // parent scan
            Assert.That(File.Exists(Path.Combine(outputFolder, @"3_SLGKVGTR(4).pdf"))); // child scan

            // clean up resources
            metadrawLogic.CleanUpResources();
            Assert.That(!metadrawLogic.FilteredListOfPsms.Any());
            Assert.That(!metadrawLogic.PsmResultFilePaths.Any());
            Assert.That(!metadrawLogic.SpectraFilePaths.Any());

            // delete output
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void MetaDraw_GlycoSearchTaskWithChildScansTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"MetaDraw_GlycoSearchTaskTest");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\leukosialin.fasta");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\sliced_glyco_hcd_ethcd.raw");

            // run task
            CommonParameters commonParameters = new CommonParameters(dissociationType: DissociationType.HCD, ms2childScanDissociationType: DissociationType.EThcD);

            Directory.CreateDirectory(outputFolder);
            var glycoSearchTask = new GlycoSearchTask() { CommonParameters = commonParameters };
            glycoSearchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, new List<string> { spectraFile }, "");

            var psmFile = Path.Combine(outputFolder, @"oglyco.psmtsv");

            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());
            Assert.That(metadrawLogic.FilteredListOfPsms.Any());

            // test results filter
            MetaDrawSettings.QValueFilter = 0.01;
            MetaDrawSettings.ShowDecoys = false;
            metadrawLogic.FilterPsms();
            Assert.That(metadrawLogic.FilteredListOfPsms.All(p => p.DecoyContamTarget == "T"));
            Assert.That(metadrawLogic.FilteredListOfPsms.All(p => p.QValue <= 0.01));

            // test text search filter (filter by full sequence)
            string filterString = @"STTAVQ";
            metadrawLogic.FilterPsmsByString(filterString);

            int c = 0;
            foreach (var filteredPsm in metadrawLogic.PeptideSpectralMatchesView)
            {
                var psmObj = (PsmFromTsv)filteredPsm;
                Assert.That(psmObj.FullSequence.Contains(filterString));
                c++;
            }
            Assert.That(c > 0);

            // test text search filter (filter by MS2 scan number)
            filterString = @"2";
            metadrawLogic.FilterPsmsByString(filterString);

            c = 0;
            foreach (var filteredPsm in metadrawLogic.PeptideSpectralMatchesView)
            {
                var psmObj = (PsmFromTsv)filteredPsm;
                Assert.That(psmObj.Ms2ScanNumber.ToString().Contains(filterString));
                c++;
            }
            Assert.That(c > 0);

            // draw PSM
            var plotView = new OxyPlot.Wpf.PlotView();
            var canvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            var psm = metadrawLogic.FilteredListOfPsms.First();

            metadrawLogic.DisplaySpectrumMatch(plotView, canvas, psm, parentChildView, out errors);
            Assert.That(errors == null || !errors.Any());

            // test that plot was drawn
            var peak = (LineSeries)plotView.Model.Series[0]; // the first m/z peak
            var peakPoints = peak.Points;
            Assert.That(Math.Round(peakPoints[0].X, 2) == 101.07); // m/z
            Assert.That(Math.Round(peakPoints[1].X, 2) == 101.07);
            Assert.That((int)peakPoints[0].Y == 0); // intensity
            Assert.That((int)peakPoints[1].Y == 3847);

            var plotAxes = plotView.Model.Axes;
            Assert.That(plotAxes.Count == 2);

            // test that base sequence annotation was drawn
            Assert.That(canvas.Children.Count > 0);

            // test that the plots were drawn in the parent/child view
            Assert.That(parentChildView.Plots.Count == 2);

            // test parent scan
            var parentPlot = parentChildView.Plots[0];
            Assert.That(parentPlot.SpectrumLabel == "Scan: 27 Dissociation Type: HCD MsOrder: 2 Selected Mz: 924.12 Retention Time: 32.65");
            int numAnnotatedResidues = psm.BaseSeq.Length;
            int numAnnotatedIons = psm.MatchedIons.Count(p => p.NeutralTheoreticalProduct.ProductType != ProductType.M
                && p.NeutralTheoreticalProduct.ProductType != ProductType.D);
            int numAnnotatedMods = psm.FullSequence.Count(p => p == '[');
            Assert.That(parentPlot.TheCanvas.Children.Count == numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods);

            peak = (LineSeries)parentPlot.Plot.Model.Series[0]; // the first m/z peak
            peakPoints = peak.Points;
            Assert.That(Math.Round(peakPoints[0].X, 2) == 101.07); // m/z
            Assert.That(Math.Round(peakPoints[1].X, 2) == 101.07);
            Assert.That((int)peakPoints[0].Y == 0); // intensity
            Assert.That((int)peakPoints[1].Y == 3847);

            // test child scan
            var childPlot = parentChildView.Plots[1];
            Assert.That(childPlot.SpectrumLabel == "Scan: 30 Dissociation Type: ETD MsOrder: 2 Selected Mz: 924.12 RetentionTime: 32.66");
            Assert.That(childPlot.TheCanvas.Children.Count > 0);
            numAnnotatedResidues = psm.BaseSeq.Length;
            numAnnotatedIons = psm.ChildScanMatchedIons[30]
                .Count(p => p.NeutralTheoreticalProduct.ProductType != ProductType.M
                && p.NeutralTheoreticalProduct.ProductType != ProductType.D);
            numAnnotatedMods = psm.FullSequence.Count(p => p == '[');
            Assert.That(childPlot.TheCanvas.Children.Count == numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods);

            peak = (LineSeries)childPlot.Plot.Model.Series[0]; // the first m/z peak
            peakPoints = peak.Points;
            Assert.That(Math.Round(peakPoints[0].X, 2) == 126.06); // m/z
            Assert.That(Math.Round(peakPoints[1].X, 2) == 126.06);
            Assert.That((int)peakPoints[0].Y == 0); // intensity
            Assert.That((int)peakPoints[1].Y == 8496);

            // write pdf
            var psmsToExport = metadrawLogic.FilteredListOfPsms.Where(p => p.FullSequence == "STTAVQTPTSGEPLVST[O-Glycosylation:H1N1 on X]SEPLSSK").ToList();
            metadrawLogic.ExportToPdf(plotView, canvas, psmsToExport, parentChildView, outputFolder, out errors);

            // test that pdf exists
            Assert.That(File.Exists(Path.Combine(outputFolder, @"27_STTAVQTPTSGEPLVST[O-Glycosylat.pdf"))); // parent scan
            Assert.That(File.Exists(Path.Combine(outputFolder, @"30_STTAVQTPTSGEPLVST[O-Glycosylat.pdf"))); // child scan

            // clean up resources
            metadrawLogic.CleanUpResources();
            Assert.That(!metadrawLogic.FilteredListOfPsms.Any());
            Assert.That(!metadrawLogic.PsmResultFilePaths.Any());
            Assert.That(!metadrawLogic.SpectraFilePaths.Any());

            // delete output
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestMetaDrawErrors()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMetaDrawErrors");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");

            Directory.CreateDirectory(outputFolder);

            // run search task
            var searchtask = new SearchTask();
            searchtask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, new List<string> { spectraFile }, "");

            var psmFile = Path.Combine(outputFolder, @"AllPSMs.psmtsv");

            // load results into metadraw (skipping spectra file, to produce an error msg)
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.PsmResultFilePaths.Add(psmFile);

            // this should produce an error because an expected spectra file is not present
            var errors = metadrawLogic.LoadFiles(loadSpectra: true, loadPsms: true);
            Assert.That(errors.Any());
            Assert.That(!metadrawLogic.FilteredListOfPsms.Any());

            // this should not produce an error because we said not to load spectra
            errors = metadrawLogic.LoadFiles(loadSpectra: false, loadPsms: true);
            Assert.That(!errors.Any());

            var psmsFromTsv = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            var plotView = new OxyPlot.Wpf.PlotView();
            var canvas = new Canvas();
            var parentChildScanPlotsView = new ParentChildScanPlotsView();

            // plotting PSM should produce an error because spectra are not loaded
            metadrawLogic.DisplaySpectrumMatch(plotView, canvas, psmsFromTsv.First(), parentChildScanPlotsView, out errors);
            Assert.That(errors.Any());

            // export to PDF should produce an error because spectra are not loaded
            metadrawLogic.ExportToPdf(plotView, canvas, new List<PsmFromTsv> { psmsFromTsv.First() }, parentChildScanPlotsView, outputFolder, out errors);
            Assert.That(errors.Any());

            // clean up resources
            metadrawLogic.CleanUpResources();

            // delete output
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestMetaDrawLoadingWithWeirdFileNames()
        {
            // test loading when the file has a periods, commas, spaces in the name
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMetaDrawLoadingWithWeirdFileNames");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");

            string pathWithPeriodInIt = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\S.m,al. lC,al.ib r.at,i ble_Ye.ast.mzML");
            File.Copy(spectraFile, pathWithPeriodInIt, true);
            spectraFile = pathWithPeriodInIt;

            Directory.CreateDirectory(outputFolder);

            // run search task
            var searchtask = new SearchTask();
            searchtask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, new List<string> { spectraFile }, "");

            var psmFile = Path.Combine(outputFolder, @"AllPSMs.psmtsv");

            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
            metadrawLogic.SpectraFilePaths.Add(pathWithPeriodInIt);
            var errors = metadrawLogic.LoadFiles(true, true);
            Assert.That(!errors.Any());

            Assert.That(metadrawLogic.FilteredListOfPsms.First().FileNameWithoutExtension == "S.m,al. lC,al.ib r.at,i ble_Ye.ast");

            var plotView = new OxyPlot.Wpf.PlotView();
            var canvas = new Canvas();
            var parentChildScanPlotsView = new ParentChildScanPlotsView();

            // plot PSM
            metadrawLogic.DisplaySpectrumMatch(plotView, canvas, metadrawLogic.FilteredListOfPsms.First(), parentChildScanPlotsView, out errors);
            Assert.That(errors == null || !errors.Any());

            // export to PDF
            metadrawLogic.ExportToPdf(plotView, canvas, new List<PsmFromTsv> { metadrawLogic.FilteredListOfPsms.First() }, parentChildScanPlotsView, outputFolder, out errors);
            Assert.That(!errors.Any());

            // clean up resources
            metadrawLogic.CleanUpResources();

            // delete output
            File.Delete(pathWithPeriodInIt);
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestMetaDrawWithSpectralLibrary()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMetaDrawWithSpectraLibrary");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch\P16858.fasta");
            string library1 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch\P16858_target.msp");
            string library2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch\P16858_decoy.msp");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch\slicedMouse.raw");

            Directory.CreateDirectory(outputFolder);

            // run search task
            var searchtask = new SearchTask();
            searchtask.RunTask(outputFolder,
                new List<DbForTask>
                {
                    new DbForTask(proteinDatabase, false),
                    new DbForTask(library1, false),
                    new DbForTask(library2, false),
                },
                new List<string> { spectraFile }, "");

            var psmFile = Path.Combine(outputFolder, @"AllPSMs.psmtsv");

            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
            metadrawLogic.SpectralLibraryPaths.Add(library1);
            metadrawLogic.SpectralLibraryPaths.Add(library2);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());

            // draw PSM
            var plotView = new OxyPlot.Wpf.PlotView();
            var canvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            var psm = metadrawLogic.FilteredListOfPsms.First();

            metadrawLogic.DisplaySpectrumMatch(plotView, canvas, psm, parentChildView, out errors);
            Assert.That(errors == null || !errors.Any());

            // test that plot was drawn
            var plotSeries = plotView.Model.Series;

            // test that library peaks were drawn in the mirror plot (these peaks have negative intensities)
            var mirrorPlotPeaks = plotSeries.Where(p => ((LineSeries)p).Points[1].Y < 0).ToList();
            Assert.That(mirrorPlotPeaks.Count == 52);

            var plotAxes = plotView.Model.Axes;
            Assert.That(plotAxes.Count == 2);

            // write pdf
            var psmsToExport = metadrawLogic.FilteredListOfPsms.Where(p => p.FullSequence == "VIHDNFGIVEGLMTTVHAITATQK").Take(1).ToList();
            metadrawLogic.ExportToPdf(plotView, canvas, psmsToExport, parentChildView, outputFolder, out errors);

            // test that pdf exists
            Assert.That(File.Exists(Path.Combine(outputFolder, @"6_VIHDNFGIVEGLMTTVHAITATQK.pdf")));

            // clean up resources
            metadrawLogic.CleanUpResources();
            Assert.That(!metadrawLogic.FilteredListOfPsms.Any());
            Assert.That(!metadrawLogic.PsmResultFilePaths.Any());
            Assert.That(!metadrawLogic.SpectraFilePaths.Any());
            Assert.That(!metadrawLogic.SpectralLibraryPaths.Any());

            // delete output
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestPsmFromTsvIonParsing()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPsmFromTsvIonParsing");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch\P16858.fasta");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch\slicedMouse.raw");

            Directory.CreateDirectory(outputFolder);

            // run search task
            var searchtask = new SearchTask();
            searchtask.RunTask(outputFolder,
                new List<DbForTask>
                {
                    new DbForTask(proteinDatabase, false),
                },
                new List<string> { spectraFile }, "");

            var psmFile = Path.Combine(outputFolder, @"AllPSMs.psmtsv");

            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());

            // assert PsmFromTsv matched ion properties
            var lines = File.ReadAllLines(psmFile);
            int ind = Array.IndexOf(lines[0].Split('\t'), "Matched Ion Mass-To-Charge Ratios");
            var ionStrings = lines[1].Split('\t')[ind].Split(new char[] { ',', ';' })
                .Select(p => p.Trim().Replace("[", string.Empty).Replace("]", string.Empty)).ToList();

            var parsedIons = metadrawLogic.FilteredListOfPsms[0].MatchedIons;

            for (int i = 0; i < ionStrings.Count; i++)
            {
                var ionString = ionStrings[i];
                var parsedIon = parsedIons[i];

                var split = ionString.Split(new char[] { '+', ':' });

                string ion = split[0];
                int charge = int.Parse(split[1]);
                double mz = double.Parse(split[2]);

                Assert.That(mz == parsedIon.Mz);
                Assert.That(mz.ToMass(charge) == parsedIon.NeutralTheoreticalProduct.NeutralMass);
                Assert.That(charge == parsedIon.Charge);
                Assert.That(ion == parsedIon.NeutralTheoreticalProduct.ProductType.ToString() + parsedIon.NeutralTheoreticalProduct.FragmentNumber);
            }

            // delete output
            Directory.Delete(outputFolder, true);
        }
    }
}