using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text.RegularExpressions;
using System.Threading;
using System.Windows;
using System.Windows.Controls;
using Chemistry;
using Easy.Common.Extensions;
using EngineLayer;
using GuiFunctions;
using GuiFunctions.MetaDraw;
using MassSpectrometry;
using NUnit.Framework;
using OxyPlot.Series;
using Omics.Fragmentation;
using OxyPlot;
using Proteomics.ProteolyticDigestion;
using Readers;
using TaskLayer;
using OxyPlot.Wpf;
using LineSeries = OxyPlot.Series.LineSeries;
using OxyPlot;

namespace Test.MetaDraw
{
    [TestFixture, Apartment(ApartmentState.STA)]
    [TestFixture]
    public class MetaDrawTest
    {
        [SetUp]
        public static void SetUp()
        {
            MetaDrawSettings.ResetSettings();
        }

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

            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);

            Assert.That(parsedPsms.Count, Is.EqualTo(10));
            Assert.That(warnings.Count, Is.EqualTo(0));

            Directory.Delete(folderPath, true);
        }


        [Test]
        public static void TestMetaDrawReadPsmFileWithoutSpectralAngle()
        {
            //test if the reader still works when psm file doesn't contain spectral angle as header.
            string noSA = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\noSAreaderTest.psmtsv");
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(noSA, out var warnings);


            Assert.That(parsedPsms.Count, Is.EqualTo(15));
            Assert.That(warnings.Count, Is.EqualTo(0));
        }

        [Test]
        public static void TestMetaDrawReadSAwhenReadingPsmFile()
        {
            SearchTask searchTask = new SearchTask();

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMetaDrawReadPsmFile");

            DbForTask db = new DbForTask(myDatabase, false);
            Directory.CreateDirectory(folderPath);

            searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "metadraw");
            string psmFile = Directory.GetFiles(folderPath).First(f => f.Contains("AllPSMs.psmtsv"));

            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);

            Assert.That(parsedPsms.First().SpectralAngle, Is.EqualTo(-1));
            Assert.That(parsedPsms.Count, Is.EqualTo(10));
            Assert.That(warnings.Count, Is.EqualTo(0));

            Directory.Delete(folderPath, true);
        }


        [Test]
        public static void TestParenthesesRemovalForSilac()
        {
            string baseSequence = "ASDF(+8.01)ASDF";
            string cleanedSequence = PsmFromTsv.RemoveParentheses(baseSequence);
            Assert.That(cleanedSequence.Equals("ASDFASDF"));
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

            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);

            Directory.Delete(folderPath, true);
        }

        [Test]
        public static void TestFindVariantCrossingIons()
        {
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\VariantCrossTest.psmtsv");
            List<string> warnings = new List<string>();
            List<PsmFromTsv> psms;

            psms = SpectrumMatchTsvReader.ReadPsmTsv(myFile, out warnings);  // test will fail if the order of psms is changed to something other than top to bottom row

            // check that variant psm properties are being parsed correctly
            Assert.That(psms[0].IdentifiedSequenceVariations, Is.EqualTo(""));
            Assert.That(psms[1].IdentifiedSequenceVariations, Is.EqualTo("A147T"));

            Assert.That(psms[0].SpliceSites, Is.EqualTo("541-541"));
            Assert.That(psms[1].SpliceSites, Is.EqualTo(""));

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
                    Assert.That(actualIons.Contains(expectedIon),
                       "VariantCrossingIons should contain ion " + expectedIon + " in file " + psms[i].FileNameWithoutExtension + ".");
                foreach (string actualIon in actualIons)
                    Assert.That(expected[i].Contains(actualIon),
                        "VariantCrossingIons should not contain ion " + actualIon + " in file " + psms[i].FileNameWithoutExtension + ".");
            }
        }

        [Test]
        public static void MetaDraw_TestStationarySequencePositioning()
        {
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.SpectralMatchResultFilePaths.Add(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\TDGPTMDSearchResults.psmtsv"));
            metadrawLogic.SpectraFilePaths.Add(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\TDGPTMDSearchSingleSpectra.mzML"));
            var errors = metadrawLogic.LoadFiles(true, true);
            Assert.That(errors.Count == 1); // Singular error should be from not loading in the rest of the spectra that the search came from

            var psm = metadrawLogic.FilteredListOfPsms.First();
            var plotView = new OxyPlot.Wpf.PlotView();
            var stationaryCanvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var sequenceAnnotationCanvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            MetaDrawSettings.FirstAAonScreenIndex = 0;
            MetaDrawSettings.NumberOfAAOnScreen = 20; // Will be dynamic based upon window size, 20 is an arbitraty number used for testing purposes
            MetaDrawSettings.DrawMatchedIons = true;
            int numAnnotatedResidues = psm.BaseSeq.Length;
            int numAnnotatedIons = psm.MatchedIons.Count;
            int numAnnotatedMods = psm.FullSequence.Count(p => p == '[');
            var peptide = new PeptideWithSetModifications(psm.FullSequence, GlobalVariables.AllModsKnownDictionary);

            // Iterates through the psm, simulating scrolling, until the sequence is scrolled as far as allowed
            for (; MetaDrawSettings.FirstAAonScreenIndex < psm.BaseSeq.Length - MetaDrawSettings.NumberOfAAOnScreen; MetaDrawSettings.FirstAAonScreenIndex++)
            {
                metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, psm);
                metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out errors);
                Assert.That(errors == null || !errors.Any());
                // Checks to see if scrollable sequence is the same each time
                Assert.That(metadrawLogic.ScrollableSequence.SequenceDrawingCanvas.Children.Count == numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods);
                Assert.That(sequenceAnnotationCanvas.Children.Count == numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods);

                // Checks to see if the stationary sequence updated with the new positioning
                string modifiedBaseSeq = psm.BaseSeq.Substring(MetaDrawSettings.FirstAAonScreenIndex, MetaDrawSettings.NumberOfAAOnScreen);
                string fullSequence = modifiedBaseSeq;
                var modDictionary = peptide.AllModsOneIsNterminus.Where(p => p.Key - 1 >= MetaDrawSettings.FirstAAonScreenIndex
                    && p.Key - 1 < (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen)).OrderByDescending(p => p.Key);
                foreach (var mod in modDictionary)
                {
                    // if modification is within the visible region
                    fullSequence = fullSequence.Insert(mod.Key - 1 - MetaDrawSettings.FirstAAonScreenIndex, "[" + mod.Value.ModificationType + ":" + mod.Value.IdWithMotif + "]");
                }

                List<MatchedFragmentIon> matchedIons = psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.ResiduePosition > MetaDrawSettings.FirstAAonScreenIndex &&
                                                       p.NeutralTheoreticalProduct.ResiduePosition < (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen)).ToList();
                int psmStartResidue = int.Parse(psm.StartAndEndResiduesInParentSequence.Split("to")[0].Replace("[", ""));
                var startAA = (MetaDrawSettings.FirstAAonScreenIndex + psmStartResidue).ToString();
                var endAA = (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen + psmStartResidue - 1).ToString();


                var expectedBaseSequence = modifiedBaseSeq.Length;
                var expectedIonAnnotations = matchedIons.Count;
                var expectedModCount = fullSequence.Count(p => p == '[');
                var expectedNumberCount = startAA.Length + endAA.Length;
                var expectedNumberLineConnectorCount = 2;

                var expected = expectedBaseSequence + expectedIonAnnotations + expectedModCount + expectedNumberCount + expectedNumberLineConnectorCount;
                Assert.That(metadrawLogic.StationarySequence.SequenceDrawingCanvas.Children.Count, Is.EqualTo(expected));
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
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
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

            MetaDrawSettings.AmbiguityFilter = "2A";
            metadrawLogic.FilterPsms();
            Assert.That(metadrawLogic.FilteredListOfPsms.All(p => p.AmbiguityLevel == "2A"));

            MetaDrawSettings.AmbiguityFilter = "No Filter";
            metadrawLogic.FilterPsms();

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

            // test text search filter (filter by organism name)
            filterString = "Sacc";
            metadrawLogic.FilterPsmsByString(filterString);

            c = 0;
            foreach (var filteredPsm in metadrawLogic.PeptideSpectralMatchesView)
            {
                var psmObj = (PsmFromTsv)filteredPsm;
                Assert.That(psmObj.OrganismName.Contains(filterString));
                c++;
            }
            Assert.That(c > 0);

            // test text search filter (filter by protein name)
            filterString = "tRNA";
            metadrawLogic.FilterPsmsByString(filterString);

            c = 0;
            foreach (var filteredPsm in metadrawLogic.PeptideSpectralMatchesView)
            {
                var psmObj = (PsmFromTsv)filteredPsm;
                Assert.That(psmObj.ProteinName.Contains(filterString));
                c++;
            }
            Assert.That(c > 0);


            // draw PSM
            var plotView = new OxyPlot.Wpf.PlotView() {Name = "plotView"};
            var stationaryCanvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var sequenceAnnotationCanvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            var psm = metadrawLogic.FilteredListOfPsms.First();

            MetaDrawSettings.FirstAAonScreenIndex = 0;
            MetaDrawSettings.NumberOfAAOnScreen = psm.BaseSeq.Length;
            MetaDrawSettings.DrawMatchedIons = true;
            MetaDrawSettings.DisplayInternalIons = false;
            metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, psm);
            metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out errors);
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

            // test with different drawing settings
            MetaDrawSettings.AnnotateCharges = true;
            MetaDrawSettings.AnnotateMzValues = true;
            metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, psm);
            metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out errors);
            Assert.That(errors == null || !errors.Any());

            MetaDrawSettings.DisplayIonAnnotations = false;
            metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, psm);
            metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out errors);
            Assert.That(errors == null || !errors.Any());

            // test that scrollable sequence annotation was drawn
            int numAnnotatedResidues = psm.BaseSeq.Length;
            int numAnnotatedIons = psm.MatchedIons.Count;
            int numAnnotatedMods = psm.FullSequence.Count(p => p == '[');
            Assert.That(metadrawLogic.ScrollableSequence.SequenceDrawingCanvas.Children.Count == numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods);
            Assert.That(sequenceAnnotationCanvas.Children.Count == numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods);

            // test that the stationary sequence annotation was drawn
            string modifiedBaseSeq = psm.BaseSeq.Substring(MetaDrawSettings.FirstAAonScreenIndex, MetaDrawSettings.NumberOfAAOnScreen);
            string fullSequence = psm.BaseSeq;
            Dictionary<int, List<string>> modDictionary = PsmFromTsv.ParseModifications(psm.FullSequence);
            foreach (var mod in modDictionary.OrderByDescending(p => p.Key))
            {
                // if modification is within the visible region
                if (mod.Key >= MetaDrawSettings.FirstAAonScreenIndex && mod.Key < (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen))
                {
                    // account for multiple modifications on the same amino acid
                    for (int i = mod.Value.Count - 1; i > -1; i--)
                    {
                        fullSequence = fullSequence.Insert(mod.Key - MetaDrawSettings.FirstAAonScreenIndex, "[" + mod.Value[i] + "]");
                        if (i >= 1)
                        {
                            fullSequence = fullSequence.Insert(mod.Key, "|");
                        }
                    }
                }
            }
            List<MatchedFragmentIon> matchedIons = psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.ResiduePosition > MetaDrawSettings.FirstAAonScreenIndex &&
                                                   p.NeutralTheoreticalProduct.ResiduePosition < (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen)).ToList();
            Assert.That(metadrawLogic.StationarySequence.SequenceDrawingCanvas.Children.Count == modifiedBaseSeq.Length + matchedIons.Count + fullSequence.Count(p => p == '[') + 
                psm.StartAndEndResiduesInParentSequence.Replace("[","").Replace("]","").Replace("to", "").Replace(" ", "").Length + 2);

            // write pdf
            var psmsToExport = metadrawLogic.FilteredListOfPsms.Where(p => p.FullSequence == "QIVHDSGR").Take(3).ToList();
            MetaDrawSettings.NumberOfAAOnScreen = psmsToExport.First().BaseSeq.Length;
            metadrawLogic.ExportPlot(plotView, stationaryCanvas, psmsToExport, parentChildView, outputFolder, out errors);

            // test that pdf exists
            Assert.That(File.Exists(Path.Combine(outputFolder, @"116_QIVHDSGR.pdf")));
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR.pdf")));
            Assert.That(File.Exists(Path.Combine(outputFolder, @"127_QIVHDSGR.pdf")));

            // test displaying a PSM with a mod
            var modPsm = metadrawLogic.FilteredListOfPsms.First(p => p.FullSequence.Contains("["));
            MetaDrawSettings.NumberOfAAOnScreen = modPsm.BaseSeq.Length;
            metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, modPsm);
            metadrawLogic.DisplaySpectrumMatch(plotView, modPsm, parentChildView, out errors);
            Assert.That(errors == null || !errors.Any());
            Assert.That(metadrawLogic.ScrollableSequence.SequenceDrawingCanvas.Children.Count == modPsm.BaseSeq.Length + modPsm.MatchedIons.Count + modPsm.FullSequence.Count(p => p == '['));
            Assert.That(sequenceAnnotationCanvas.Children.Count == modPsm.BaseSeq.Length + modPsm.MatchedIons.Count + modPsm.FullSequence.Count(p => p == '['));
            modifiedBaseSeq = modPsm.BaseSeq.Substring(MetaDrawSettings.FirstAAonScreenIndex, MetaDrawSettings.NumberOfAAOnScreen);
            fullSequence = modPsm.BaseSeq;
            modDictionary = PsmFromTsv.ParseModifications(modPsm.FullSequence);
            foreach (var mod in modDictionary.OrderByDescending(p => p.Key))
            {
                // if modification is within the visible region
                if (mod.Key >= MetaDrawSettings.FirstAAonScreenIndex && mod.Key < (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen))
                {
                    // account for multiple modifications on the same amino acid
                    for (int i = mod.Value.Count - 1; i > -1; i--)
                    {
                        fullSequence = fullSequence.Insert(mod.Key - MetaDrawSettings.FirstAAonScreenIndex, "[" + mod.Value[i] + "]");
                        if (i >= 1)
                        {
                            fullSequence = fullSequence.Insert(mod.Key, "|");
                        }
                    }
                }
            }
            matchedIons = modPsm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.ResiduePosition > MetaDrawSettings.FirstAAonScreenIndex &&
                                                p.NeutralTheoreticalProduct.ResiduePosition < (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen)).ToList();
            Assert.That(metadrawLogic.StationarySequence.SequenceDrawingCanvas.Children.Count == modifiedBaseSeq.Length + matchedIons.Count + fullSequence.Count(p => p == '[') + 
                psm.StartAndEndResiduesInParentSequence.Replace("[", "").Replace("]", "").Replace("to", "").Replace(" ", "").Length + 2);


            // get scan from psm
            var scan = metadrawLogic.GetMs2ScanFromPsm(psm);
            Assert.That(scan.OneBasedScanNumber, Is.EqualTo(psm.Ms2ScanNumber));
            Assert.That(scan.MsnOrder, Is.EqualTo(2));

            // clean up resources
            metadrawLogic.CleanUpResources();
            Assert.That(!metadrawLogic.FilteredListOfPsms.Any());
            Assert.That(!metadrawLogic.SpectralMatchResultFilePaths.Any());
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
            metadrawLogic.SpectralMatchResultFilePaths.Add(csmFile);
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
            var plotView = new OxyPlot.Wpf.PlotView() { Name = "plotView" };
            var stationaryCanvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var sequenceAnnotationCanvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            var csm = metadrawLogic.FilteredListOfPsms.First();

            metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, csm);
            metadrawLogic.DisplaySpectrumMatch(plotView, csm, parentChildView, out errors);
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
            Assert.That(stationaryCanvas.Children.Count > 0);

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
            numAnnotatedIons = csm.ChildScanMatchedIons[3].Concat((csm as PsmFromTsv).BetaPeptideChildScanMatchedIons[3])
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
            metadrawLogic.ExportPlot(plotView, metadrawLogic.StationarySequence.SequenceDrawingCanvas, psmsToExport, parentChildView, outputFolder, out errors);

            // write pdf with legend
            Canvas ptmLegend = new();
            System.Windows.Size legendSize = new(100, 100);
            ptmLegend.Measure(legendSize);
            ptmLegend.Arrange(new Rect(legendSize));
            ptmLegend.UpdateLayout();
            Vector ptmLegendVector = new(10, 10);
            metadrawLogic.ExportPlot(plotView, metadrawLogic.StationarySequence.SequenceDrawingCanvas, psmsToExport,
                parentChildView, outputFolder, out errors, ptmLegend, ptmLegendVector);

            // test that pdf exists
            Assert.That(File.Exists(Path.Combine(outputFolder, @"2_SLGKVGTR(4)EKVLTSSAR(2).pdf"))); // parent scan
            Assert.That(File.Exists(Path.Combine(outputFolder, @"3_SLGKVGTR(4)EKVLTSSAR(2).pdf"))); // child scan

            // clean up resources
            metadrawLogic.CleanUpResources();
            Assert.That(!metadrawLogic.FilteredListOfPsms.Any());
            Assert.That(!metadrawLogic.SpectralMatchResultFilePaths.Any());
            Assert.That(!metadrawLogic.SpectraFilePaths.Any());

            // delete output
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestMetaDrawXlSpectralLibrary()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlOutputTestFile");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\RibosomeGO.fasta");
            string library1 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SpectralLibrarySearch\CrosslinkSpectralLibrary.msp");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\2017-11-21_XL_DSSO_Ribosome_RT60min_28800-28898.mzML");

            Directory.CreateDirectory(outputFolder);

            XLSearchTask xLSearch = new XLSearchTask
            {
                XlSearchParameters = new XlSearchParameters
                {
                    WriteSpectralLibrary = true,
                    CrosslinkAtCleavageSite = true
                }
            };
            xLSearch.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) },
                new List<string> { spectraFile }, "test");

            var psmFile = Path.Combine(outputFolder, @"XL_Intralinks.tsv");

            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
            metadrawLogic.SpectralLibraryPaths.Add(library1);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());

            // draw PSM
            var plotView = new OxyPlot.Wpf.PlotView() { Name = "plotView" };
            var canvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var stationaryCanvas = new Canvas();
            var sequenceAnnotationCanvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            var psm = metadrawLogic.FilteredListOfPsms.First();

            MsDataFile file = new Mzml(spectraFile);

            MetaDrawSettings.FirstAAonScreenIndex = 0;
            MetaDrawSettings.NumberOfAAOnScreen = metadrawLogic.FilteredListOfPsms.First().BaseSeq.Length;
            metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, psm);
            metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out errors);
            Assert.That(errors == null || !errors.Any());

            // test that plot was drawn
            var plotSeries = plotView.Model.Series;

            // test that library peaks were drawn in the mirror plot (these peaks have negative intensities)
            var mirrorPlotPeaks = plotSeries.Where(p => ((LineSeries)p).Points[1].Y < 0).ToList();
            Assert.That(mirrorPlotPeaks.Count, Is.EqualTo(59));

            var plotAxes = plotView.Model.Axes;
            Assert.That(plotAxes.Count == 2);

            // write pdf
            var psmsToExport = metadrawLogic.FilteredListOfPsms.Where(p => (p as PsmFromTsv).UniqueSequence == "LLDNAAADLAAISGQKPLITKAR(21)ITLNMGVGEAIADKK(14)").Take(1).ToList();
            metadrawLogic.ExportPlot(plotView, canvas, psmsToExport, parentChildView, outputFolder, out errors);


            // write pdf with legend
            Canvas ptmLegend = new();
            System.Windows.Size legendSize = new(100, 100);
            ptmLegend.Measure(legendSize);
            ptmLegend.Arrange(new Rect(legendSize));
            ptmLegend.UpdateLayout();
            Vector ptmLegendVector = new(10, 10);
            metadrawLogic.ExportPlot(plotView, metadrawLogic.StationarySequence.SequenceDrawingCanvas, psmsToExport,
                parentChildView, outputFolder, out errors, ptmLegend, ptmLegendVector);

            // test that pdf exists
            Assert.That(File.Exists(Path.Combine(outputFolder, @"13_LLDNAAADLAAISGQKPLITKAR(21)ITL.pdf"))); // Name can only be 30  characters long

            // clean up resources
            metadrawLogic.CleanUpResources();
            Assert.That(!metadrawLogic.FilteredListOfPsms.Any());
            Assert.That(!metadrawLogic.SpectralMatchResultFilePaths.Any());
            Assert.That(!metadrawLogic.SpectraFilePaths.Any());
            Assert.That(!metadrawLogic.SpectralLibraryPaths.Any());

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
            glycoSearchTask._glycoSearchParameters.OxoniumIonFilt = false; //turn off the diagnostic filter, because the case we use have 272,294 oxonium ions but assigned the N1H1 to that.
            glycoSearchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, new List<string> { spectraFile }, "");

            var psmFile = Path.Combine(outputFolder, @"oglyco.psmtsv");

            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
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
            var plotView = new OxyPlot.Wpf.PlotView(){Name = "plotView"};
            var parentChildView = new ParentChildScanPlotsView();
            var stationaryCanvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var sequenceAnnotationCanvas = new Canvas();
            var psm = metadrawLogic.FilteredListOfPsms.First();
            MetaDrawSettings.NumberOfAAOnScreen = psm.BaseSeq.Length;

            metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, psm);
            metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out errors);
            Assert.That(errors == null || !errors.Any());

            // test that plot was drawn
            var peak = (LineSeries)plotView.Model.Series[0]; // the first m/z peak
            var peakPoints = peak.Points;
            Assert.That(Math.Round(peakPoints[0].X, 2) == 101.07); // m/z
            Assert.That(Math.Round(peakPoints[1].X, 2) == 101.07);
            Assert.That((int)peakPoints[0].Y == 0); // intensity
            Assert.That((int)peakPoints[1].Y == 3847);

            var plotAxes = plotView.Model.Axes;
            Assert.That(plotAxes.Count, Is.EqualTo(2));

            // test that base sequence annotation was drawn
            Assert.That(metadrawLogic.StationarySequence.SequenceDrawingCanvas.Children.Count > 0);

            // test that the plots were drawn in the parent/child view
            Assert.That(parentChildView.Plots.Count, Is.EqualTo(2));

            // test parent scan
            var parentPlot = parentChildView.Plots[0];
            Assert.That(parentPlot.SpectrumLabel, Is.EqualTo("Scan: 27 Dissociation Type: HCD MsOrder: 2 Selected Mz: 924.12 Retention Time: 32.65"));
            int numAnnotatedResidues = psm.BaseSeq.Length;
            int numAnnotatedIons = psm.MatchedIons.Count(p => p.NeutralTheoreticalProduct.ProductType != ProductType.M
                && p.NeutralTheoreticalProduct.ProductType != ProductType.D);
            int numAnnotatedMods = psm.FullSequence.Count(p => p == '[');
            Assert.That(parentPlot.TheCanvas.Children.Count, Is.EqualTo(numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods));

            peak = (LineSeries)parentPlot.Plot.Model.Series[0]; // the first m/z peak
            peakPoints = peak.Points;
            Assert.That(Math.Round(peakPoints[0].X, 2), Is.EqualTo(101.07)); // m/z
            Assert.That(Math.Round(peakPoints[1].X, 2), Is.EqualTo(101.07));
            Assert.That((int)peakPoints[0].Y, Is.EqualTo(0)); // intensity
            Assert.That((int)peakPoints[1].Y, Is.EqualTo(3847));

            // test child scan
            var childPlot = parentChildView.Plots[1];
            Assert.That(childPlot.SpectrumLabel, Is.EqualTo("Scan: 30 Dissociation Type: EThcD MsOrder: 2 Selected Mz: 924.12 RetentionTime: 32.66"));
            Assert.That(childPlot.TheCanvas.Children.Count > 0);
            numAnnotatedResidues = psm.BaseSeq.Length;
            numAnnotatedIons = psm.ChildScanMatchedIons[30]
                .Count(p => p.NeutralTheoreticalProduct.ProductType != ProductType.M
                && p.NeutralTheoreticalProduct.ProductType != ProductType.D);
            numAnnotatedMods = psm.FullSequence.Count(p => p == '[');
            Assert.That(childPlot.TheCanvas.Children.Count, Is.EqualTo(numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods));

            peak = (LineSeries)childPlot.Plot.Model.Series[0]; // the first m/z peak
            peakPoints = peak.Points;
            Assert.That(Math.Round(peakPoints[0].X, 2), Is.EqualTo(126.05)); // m/z
            Assert.That(Math.Round(peakPoints[1].X, 2), Is.EqualTo(126.06));
            Assert.That((int)peakPoints[0].Y, Is.EqualTo(0)); // intensity
            Assert.That((int)peakPoints[1].Y, Is.EqualTo(8496));

            // write pdf
            var psmsToExport = metadrawLogic.FilteredListOfPsms.Where(p => p.FullSequence == "STTAVQTPTSGEPLVST[O-linked glycosylation:H1N1 on T]SEPLSSK").ToList();
            metadrawLogic.ExportPlot(plotView, metadrawLogic.StationarySequence.SequenceDrawingCanvas, psmsToExport, parentChildView, outputFolder, out errors);

            // test that pdf exists
            Assert.That(File.Exists(Path.Combine(outputFolder, @"27_STTAVQTPTSGEPLVST[O-linked gly.pdf"))); // parent scan
            Assert.That(File.Exists(Path.Combine(outputFolder, @"30_STTAVQTPTSGEPLVST[O-linked gly.pdf"))); // child scan
            Directory.Delete(outputFolder, true);

            Canvas ptmLegend = new();
            System.Windows.Size legendSize = new(100, 100);
            ptmLegend.Measure(legendSize);
            ptmLegend.Arrange(new Rect(legendSize));
            ptmLegend.UpdateLayout();   
            Vector ptmLegendVector = new(10, 10);
            metadrawLogic.ExportPlot(plotView, metadrawLogic.StationarySequence.SequenceDrawingCanvas, psmsToExport, parentChildView, outputFolder, out errors, ptmLegend, ptmLegendVector);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"27_STTAVQTPTSGEPLVST[O-linked gly.pdf"))); // parent scan
            Assert.That(File.Exists(Path.Combine(outputFolder, @"30_STTAVQTPTSGEPLVST[O-linked gly.pdf"))); // child scan

            // test that a directory is created if it does not exist
            Directory.Delete(outputFolder, true);
            metadrawLogic.ExportPlot(plotView, metadrawLogic.StationarySequence.SequenceDrawingCanvas, psmsToExport, parentChildView, outputFolder, out errors);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"27_STTAVQTPTSGEPLVST[O-linked gly.pdf"))); // parent scan
            Assert.That(File.Exists(Path.Combine(outputFolder, @"30_STTAVQTPTSGEPLVST[O-linked gly.pdf"))); // child scan



            // clean up resources
            metadrawLogic.CleanUpResources();
            Assert.That(!metadrawLogic.FilteredListOfPsms.Any());
            Assert.That(!metadrawLogic.SpectralMatchResultFilePaths.Any());
            Assert.That(!metadrawLogic.SpectraFilePaths.Any());

            // delete output
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void MetaDraw_TestChimeraScanSpectrumDisplay()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"MetaDraw_SearchTaskTest");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\hela_snip_for_unitTest.fasta");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            Regex illegalInFileName = new Regex(@"[\\/:*?""<>|]");
            Directory.CreateDirectory(outputFolder);

            // run search task
            var searchtask = new SearchTask()
            {
                SearchParameters = new SearchParameters() { MinAllowedInternalFragmentLength = 4},
            };
            searchtask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, new List<string> { spectraFile }, "");

            var psmFile = Path.Combine(outputFolder, @"AllPSMs.psmtsv");

            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(true, true);
            Assert.That(!errors.Any());
            Assert.That(metadrawLogic.FilteredListOfPsms.Any());

            // fix the scan number due to the trimmed spectra file
            foreach (var psm in metadrawLogic.FilteredListOfPsms)
            {
                var type = typeof(SpectrumMatchFromTsv);
                var field = type.GetField("<Ms2ScanNumber>k__BackingField", BindingFlags.Instance | BindingFlags.NonPublic);
                field.SetValue(psm, psm.Ms2ScanNumber + 27300);

                field = type.GetField("<PrecursorScanNum>k__BackingField", BindingFlags.Instance | BindingFlags.NonPublic);
                field.SetValue(psm, psm.PrecursorScanNum + 27300);
            }

            var metaDrawDynamicScanConnection = metadrawLogic.MsDataFiles;

            var chimeraAnalysisTab = new ChimeraAnalysisTabViewModel(
                metadrawLogic.FilteredListOfPsms.ToList(), metadrawLogic.MsDataFiles,
                outputFolder);
            // test plotting on each instance of chimeras in this dataset
            var plotView = new OxyPlot.Wpf.PlotView() { Name = "chimeraPlot" };
            foreach (var chimeraGroup in chimeraAnalysisTab.ChimeraGroupViewModels)
            {
                Assert.That(chimeraGroup.Count, Is.GreaterThanOrEqualTo(2));
                MsDataScan chimericScan = metaDrawDynamicScanConnection.First().Value
                    .GetOneBasedScanFromDynamicConnection(chimeraGroup.Ms2ScanNumber);

                // plot the first chimera and test the results
                var plot = new ChimeraSpectrumMatchPlot(plotView, chimeraGroup);
                Assert.That(errors == null || !errors.Any());

                // test plot was drawn
                var model = plotView.Model;
                Assert.That(model, Is.Not.Null);
                Assert.That(plotView.Model.Axes.Count == 2);

                var peakPoints = ((LineSeries)model.Series[0]).Points;

                Assert.That(Math.Round(peakPoints[0].X, 2) == Math.Round(chimericScan.MassSpectrum.XArray[0], 2)); // m/z
                Assert.That(Math.Round(peakPoints[1].X, 2) == Math.Round(chimericScan.MassSpectrum.XArray[0], 2));
                Assert.That((int)peakPoints[0].Y == 0); // intensity
                Assert.That((int)peakPoints[1].Y == (int)chimericScan.MassSpectrum.YArray[0]);

                // Classify All Fragments
                Dictionary<OxyColor, List<MatchedFragmentIon>> ionDict = new();
                var accessionToSharedProteinColor = chimeraGroup.ChimericPsms
                    .GroupBy(p => p.Psm.Accession)
                    .Select((p, i) =>
                    (p.Key, Color: ChimeraGroupViewModel.ColorByProteinDictionary[i][0])).ToDictionary(p => p.Key, p => p.Color);
                var fullSequenceToUniqueColorDict = chimeraGroup.ChimericPsms
                    .GroupBy(p => p.Psm.Accession)
                    .SelectMany((p, i) =>
                    {
                        if (p.Count() == 1)
                        {
                            return [(p.First().Psm.FullSequence, Color: ChimeraGroupViewModel.ColorByProteinDictionary[i][1])];
                        }
                        else
                        {
                            return p.Select((m, j) => (m.Psm.FullSequence, Color: ChimeraGroupViewModel.ColorByProteinDictionary[i][j+1]));
                        }
                    }).ToDictionary(p => p.FullSequence, p => p.Color);

                foreach (var ionGroup in chimeraGroup.ChimericPsms
                    .SelectMany((p, i) => p.Psm.MatchedIons
                        .Select(ion => (ion, i, p.Psm.Accession, p.Psm.FullSequence)))
                    .GroupBy(p => p.ion.Mz))
                {
                    if (ionGroup.Count() == 1) // unshared peak
                    {
                        var color = fullSequenceToUniqueColorDict[ionGroup.First().FullSequence];
                        if (ionDict.ContainsKey(color))
                        {
                            ionDict[color].Add(ionGroup.First().ion);
                        }
                        else
                        {
                            ionDict[color] = new List<MatchedFragmentIon> { ionGroup.First().ion };
                        }
                    }
                    else // shared peak
                    {
                        int accessionCount = ionGroup.Select(p => p.Accession).Distinct().Count();
                        int fullSequenceCount = ionGroup.Select(p => p.FullSequence).Distinct().Count();
                        if (accessionCount == 1 && fullSequenceCount > 1)
                        {
                            // shared peak from different proteoforms of the same protein
                            var color = accessionToSharedProteinColor[ionGroup.First().Accession];
                            if (ionDict.ContainsKey(color))
                            {
                                ionDict[color].AddRange(ionGroup.Select(p => p.ion));
                            }
                            else
                            {
                                ionDict[color] = ionGroup.Select(p => p.ion).ToList();
                            }
                        }
                        else if (accessionCount == 1 && fullSequenceCount == 1) // shared peak from the same protein
                        {
                            var color = fullSequenceToUniqueColorDict[ionGroup.First().FullSequence];
                            if (ionDict.ContainsKey(color))
                            {
                                ionDict[color].AddRange(ionGroup.Select(p => p.ion));
                            }
                            else
                            {
                                ionDict[color] = ionGroup.Select(p => p.ion).ToList();
                            }
                        }
                        else // shared peak from different proteins
                        {
                            foreach (var ion in ionGroup)
                            {
                                if (ionDict.ContainsKey(ChimeraGroupViewModel.MultipleProteinSharedColor))
                                {
                                    ionDict[ChimeraGroupViewModel.MultipleProteinSharedColor].Add(ion.ion);
                                }
                                else
                                {
                                    ionDict[ChimeraGroupViewModel.MultipleProteinSharedColor] = [ion.ion];
                                }
                            }
                        }

                    }
                }

                // shared matched ions are default color
                var expectedSharedIons = ionDict[ChimeraGroupViewModel.MultipleProteinSharedColor].DistinctBy(p => p.Mz).ToList();
                int drawnIonsShared = model.Series.Count(p => ((LineSeries)p).Color == ChimeraGroupViewModel.MultipleProteinSharedColor);
                Assert.That(drawnIonsShared, Is.EqualTo(expectedSharedIons.Count));
                if (expectedSharedIons.Count > 0 ) 
                {
                    Assert.That(chimeraGroup.MatchedFragmentIonsByColor.ContainsKey(ChimeraGroupViewModel.MultipleProteinSharedColor));
                    Assert.That(chimeraGroup.MatchedFragmentIonsByColor[ChimeraGroupViewModel.MultipleProteinSharedColor].Count, Is.EqualTo(expectedSharedIons.Count));
                }
                else
                    Assert.That(!chimeraGroup.MatchedFragmentIonsByColor.ContainsKey(ChimeraGroupViewModel.MultipleProteinSharedColor));

                // Ions Shared by multiple proteoforms are annotated correctly
                foreach (var color in chimeraGroup.MatchedFragmentIonsByColor.Keys)
                {
                    if (color == ChimeraGroupViewModel.MultipleProteinSharedColor) continue;
                    Assert.That(ionDict.ContainsKey(color));
                    Assert.That(chimeraGroup.MatchedFragmentIonsByColor[color].Count, Is.EqualTo(ionDict[color].Count));
                }

                // test with different drawing settings
                MetaDrawSettings.AnnotateCharges = true;
                MetaDrawSettings.AnnotateMzValues = true;
                plot = new ChimeraSpectrumMatchPlot(plotView, chimeraGroup);
                Assert.That(errors == null || !errors.Any());

                MetaDrawSettings.DisplayIonAnnotations = false;
                plot = new ChimeraSpectrumMatchPlot(plotView, chimeraGroup);
                Assert.That(errors == null || !errors.Any());

                MetaDrawSettings.DisplayInternalIons = false;
                plot = new ChimeraSpectrumMatchPlot(plotView, chimeraGroup);
                Assert.That(errors == null || !errors.Any());
                MetaDrawSettings.DisplayInternalIons = true;
            }

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
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);

            // this should produce an error because an expected spectra file is not present
            var errors = metadrawLogic.LoadFiles(loadSpectra: true, loadPsms: true);
            Assert.That(errors.Any());
            Assert.That(!metadrawLogic.FilteredListOfPsms.Any());

            // this should not produce an error because we said not to load spectra
            errors = metadrawLogic.LoadFiles(loadSpectra: false, loadPsms: true);
            Assert.That(!errors.Any());

            var psmsFromTsv = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            var plotView = new OxyPlot.Wpf.PlotView();
            var canvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var parentChildScanPlotsView = new ParentChildScanPlotsView();

            // plotting PSM should produce an error because spectra are not loaded
            metadrawLogic.DisplaySpectrumMatch(plotView, psmsFromTsv.First(), parentChildScanPlotsView, out errors);
            Assert.That(errors.Any());

            // export to PDF should produce an error because spectra are not loaded
            MetaDrawSettings.NumberOfAAOnScreen = psmsFromTsv.First().BaseSeq.Length - 1;
            metadrawLogic.ExportPlot(plotView, canvas, new List<SpectrumMatchFromTsv> { psmsFromTsv.First() }, parentChildScanPlotsView, outputFolder, out errors);
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
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\S.m,al. lC,al.ib r.at,i ble_Ye.ast.mzML");


            Directory.CreateDirectory(outputFolder);

            // run search task
            var searchtask = new SearchTask();
            searchtask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, new List<string> { spectraFile }, "");

            var psmFile = Path.Combine(outputFolder, @"AllPSMs.psmtsv");

            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            var errors = metadrawLogic.LoadFiles(true, true);
            Assert.That(!errors.Any());

            Assert.That(metadrawLogic.FilteredListOfPsms.First().FileNameWithoutExtension == "S.m,al. lC,al.ib r.at,i ble_Ye.ast");

            var plotView = new OxyPlot.Wpf.PlotView();
            var canvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var stationaryCanvas = new Canvas();
            var sequenceAnnotationCanvas = new Canvas();
            var parentChildScanPlotsView = new ParentChildScanPlotsView();

            // plot PSM
            MetaDrawSettings.FirstAAonScreenIndex = 0;
            MetaDrawSettings.NumberOfAAOnScreen = metadrawLogic.FilteredListOfPsms.First().BaseSeq.Length;
            metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, metadrawLogic.FilteredListOfPsms.First());
            metadrawLogic.DisplaySpectrumMatch(plotView, metadrawLogic.FilteredListOfPsms.First(), parentChildScanPlotsView, out errors);
            Assert.That(errors == null || !errors.Any());

            // export to PDF
            metadrawLogic.ExportPlot(plotView, canvas, new List<SpectrumMatchFromTsv> { metadrawLogic.FilteredListOfPsms.First() }, parentChildScanPlotsView, outputFolder, out errors);
            Assert.That(!errors.Any());

            // clean up resources
            metadrawLogic.CleanUpResources();
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
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
            metadrawLogic.SpectralLibraryPaths.Add(library1);
            metadrawLogic.SpectralLibraryPaths.Add(library2);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());

            // draw PSM
            var plotView = new OxyPlot.Wpf.PlotView() { Name = "plotView" };
            var canvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var stationaryCanvas = new Canvas();
            var sequenceAnnotationCanvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            var psm = metadrawLogic.FilteredListOfPsms.First();

            MetaDrawSettings.FirstAAonScreenIndex = 0;
            MetaDrawSettings.NumberOfAAOnScreen = metadrawLogic.FilteredListOfPsms.First().BaseSeq.Length;
            metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, psm);
            metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out errors);
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
            metadrawLogic.ExportPlot(plotView, canvas, psmsToExport, parentChildView, outputFolder, out errors);

            // test that pdf exists
            Assert.That(File.Exists(Path.Combine(outputFolder, @"6_VIHDNFGIVEGLMTTVHAITATQK.pdf")));

            // clean up resources
            metadrawLogic.CleanUpResources();
            Assert.That(!metadrawLogic.FilteredListOfPsms.Any());
            Assert.That(!metadrawLogic.SpectralMatchResultFilePaths.Any());
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
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
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
                Assert.That(mz.ToMass(charge) - parsedIon.MassErrorDa == parsedIon.NeutralTheoreticalProduct.NeutralMass);
                Assert.That(charge == parsedIon.Charge);
                Assert.That(ion == parsedIon.NeutralTheoreticalProduct.ProductType.ToString() + parsedIon.NeutralTheoreticalProduct.FragmentNumber);
            }

            //check that the proteoform classification was correct
            Assert.That(metadrawLogic.FilteredListOfPsms[0].AmbiguityLevel.Equals("1"));

            // delete output
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void SequenceCoverageMapTest()
        {
            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\InternalTest.mgf");
            string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SequenceCoverageTestPSM.psmtsv");
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());

            // draw sequence coverage for PSM
            var sequenceText = new Canvas();
            var map = new Canvas();
            var psm = metadrawLogic.FilteredListOfPsms.First();

            metadrawLogic.DrawSequenceCoverageMap(psm, sequenceText, map);

            //test no errors
            Assert.That(errors == null || !errors.Any());
        }

        [Test]
        public static void IonsWithNoTerminusTest() //if internal fragments are selected using TerminusType == None, then "M" and "D" ions get lumped in and cause a crash
        {
            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\ms2mixed_bsa_xlink.mzML");
            string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\XL_Intralinks_MIons.tsv");
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());

            // draw sequence coverage for PSM
            var sequenceText = new Canvas();
            var map = new Canvas();
            var psm = metadrawLogic.FilteredListOfPsms.First();

            metadrawLogic.DrawSequenceCoverageMap(psm, sequenceText, map);

            //test no errors
            Assert.That(errors == null || !errors.Any());
        }

        [Test]
        public static void TestMetaDrawLogicCleanUp()
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
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
            metadrawLogic.SpectralLibraryPaths.Add(library1);
            metadrawLogic.SpectralLibraryPaths.Add(library2);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());

            // draw PSM
            var plotView = new OxyPlot.Wpf.PlotView() { Name = "plotView" };
            var canvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var stationaryCanvas = new Canvas();
            var sequenceAnnotationCanvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            var psm = metadrawLogic.FilteredListOfPsms.First();

            MetaDrawSettings.FirstAAonScreenIndex = 0;
            MetaDrawSettings.NumberOfAAOnScreen = metadrawLogic.FilteredListOfPsms.First().BaseSeq.Length;
            metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, psm);
            metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out errors);
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
            metadrawLogic.ExportPlot(plotView, canvas, psmsToExport, parentChildView, outputFolder, out errors);

            // test that pdf exists
            Assert.That(File.Exists(Path.Combine(outputFolder, @"6_VIHDNFGIVEGLMTTVHAITATQK.pdf")));

            // clean up resources
            metadrawLogic.CleanUpSpectraFiles();
            Assert.That(!metadrawLogic.SpectraFilePaths.Any());

            metadrawLogic.CleanUpPSMFiles();
            Assert.That(!metadrawLogic.FilteredListOfPsms.Any());
            Assert.That(!metadrawLogic.SpectralMatchResultFilePaths.Any());

            metadrawLogic.CleanUpSpectralLibraryFiles();
            Assert.That(!metadrawLogic.SpectralLibraryPaths.Any());

            // delete output
            Directory.Delete(outputFolder, true);

        }

        [Test]
        public static void TestMetaDrawOutputFormats()
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
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());
            Assert.That(metadrawLogic.FilteredListOfPsms.Any());

            var plotView = new OxyPlot.Wpf.PlotView() { Name = "plotView" };
            var stationaryCanvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var sequenceAnnotationCanvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            var psm = metadrawLogic.FilteredListOfPsms.Where(p => p.FullSequence == "QIVHDSGR").First();

            // drawing the first psm
            MetaDrawSettings.FirstAAonScreenIndex = 0;
            MetaDrawSettings.NumberOfAAOnScreen = psm.BaseSeq.Length;
            MetaDrawSettings.DrawMatchedIons = true;
            metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, psm);
            metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out errors);
            Assert.That(errors == null || !errors.Any());

            // export each file type and ensure they exist
            var psmsToExport = new List<SpectrumMatchFromTsv> { metadrawLogic.FilteredListOfPsms.Where(p => p.FullSequence == "QIVHDSGR").First() };
            MetaDrawSettings.NumberOfAAOnScreen = psmsToExport.First().BaseSeq.Length;
            metadrawLogic.ExportPlot(plotView, stationaryCanvas, psmsToExport, parentChildView, outputFolder, out errors);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR.pdf")));

            MetaDrawSettings.ExportType = "Png";
            metadrawLogic.ExportPlot(plotView, stationaryCanvas, psmsToExport, parentChildView, outputFolder, out errors);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR.png")));

            MetaDrawSettings.ExportType = "Jpeg";
            metadrawLogic.ExportPlot(plotView, stationaryCanvas, psmsToExport, parentChildView, outputFolder, out errors);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR.jpeg")));

            MetaDrawSettings.ExportType = "Tiff";
            metadrawLogic.ExportPlot(plotView, stationaryCanvas, psmsToExport, parentChildView, outputFolder, out errors);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR.tiff")));

            MetaDrawSettings.ExportType = "Wmf";
            metadrawLogic.ExportPlot(plotView, stationaryCanvas, psmsToExport, parentChildView, outputFolder, out errors);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR.wmf")));

            MetaDrawSettings.ExportType = "Bmp";
            metadrawLogic.ExportPlot(plotView, stationaryCanvas, psmsToExport, parentChildView, outputFolder, out errors);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR.bmp")));

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestMetaDrawCombineBitmapTryCatch()
        {
            List<Bitmap> bmList = new();
            List<System.Windows.Point> pts = new();
            Assert.Throws<ArgumentException>(() => MetaDrawLogic.CombineBitmap(bmList, pts, false));

            Bitmap bigBmp = new(500, 500, System.Drawing.Imaging.PixelFormat.Format32bppArgb);
            for (int i = 0; i < 10; i++)
            {
                bmList.Add(bigBmp);
            }
            Assert.Throws<ArgumentOutOfRangeException>(() => MetaDrawLogic.CombineBitmap(bmList, pts, false));

        }

        [Test]
        public static void EnsureNoCrashesWithNGlyco()
        {
            string psmsPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "nglyco_f5.psmtsv");
            var psmOfInterest = SpectrumMatchTsvReader.ReadPsmTsv(psmsPath, out var warnings).First(p => p.FullSequence.Contains("LLSTEGSQ"));
            Assert.That(!warnings.Any());

            MetaDrawLogic metaDrawLogic = new();
            metaDrawLogic.SpectralMatchResultFilePaths.Add(psmsPath);
            var errors = metaDrawLogic.LoadFiles(false, true);
            Assert.That(!errors.Any());
            Assert.That(metaDrawLogic.FilteredListOfPsms.Any());
            
            var stationaryCanvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var sequenceAnnotationCanvas = new Canvas();

            // drawing the first psm
            MetaDrawSettings.FirstAAonScreenIndex = 0;
            MetaDrawSettings.NumberOfAAOnScreen = psmOfInterest.BaseSeq.Length;
            MetaDrawSettings.DrawMatchedIons = true;
            metaDrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, psmOfInterest);
            Assert.That(errors == null || !errors.Any());
        }

        [Test]
        public static void TestMetaDrawSequenceDisplayOutputs()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"MetaDraw_SearchTaskTest");
            string proteinDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");

            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);
            Directory.CreateDirectory(outputFolder);

            // run search task
            var searchtask = new SearchTask();
            searchtask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(proteinDatabase, false) }, new List<string> { spectraFile }, "");

            var psmFile = Path.Combine(outputFolder, @"AllPSMs.psmtsv");

            // load results into metadraw
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.SpectraFilePaths.Add(spectraFile);
            metadrawLogic.SpectralMatchResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(true, true);

            Assert.That(!errors.Any());
            Assert.That(metadrawLogic.FilteredListOfPsms.Any());

            var plotView = new OxyPlot.Wpf.PlotView();
            var stationaryCanvas = new Canvas();
            var scrollableCanvas = new Canvas();
            var mapCanvas = new Canvas() { Width = 200, Height = 200};
            var textCanvas = new Canvas() { Width = 200, Height = 200 };
            var sequenceAnnotationCanvas = new Canvas();
            var parentChildView = new ParentChildScanPlotsView();
            var ptmLegend = new System.Windows.UIElement();
            var psm = metadrawLogic.FilteredListOfPsms.Where(p => p.FullSequence == "QIVHDSGR").First();

            // drawing the first psm
            MetaDrawSettings.FirstAAonScreenIndex = 0;
            MetaDrawSettings.NumberOfAAOnScreen = psm.BaseSeq.Length;
            MetaDrawSettings.DrawMatchedIons = true;
            metadrawLogic.DisplaySequences(stationaryCanvas, scrollableCanvas, sequenceAnnotationCanvas, psm);
            metadrawLogic.DisplaySpectrumMatch(plotView, psm, parentChildView, out errors);
            Assert.That(errors == null || !errors.Any());

            // ensure each file type can be outputted by each method
            var psmToExport = metadrawLogic.FilteredListOfPsms.Where(p => p.FullSequence == "QIVHDSGR").First();
            MetaDrawSettings.NumberOfAAOnScreen = psmToExport.BaseSeq.Length;
            metadrawLogic.ExportSequenceCoverage(textCanvas, mapCanvas, outputFolder, psmToExport);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR_SequenceCoverage.pdf")));
            metadrawLogic.ExportAnnotatedSequence(sequenceAnnotationCanvas, ptmLegend, psmToExport, outputFolder, 200);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR_SequenceAnnotation.pdf")));

            MetaDrawSettings.ExportType = "Png";
            metadrawLogic.ExportSequenceCoverage(textCanvas, mapCanvas, outputFolder, psmToExport);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR_SequenceCoverage.png")));
            metadrawLogic.ExportAnnotatedSequence(sequenceAnnotationCanvas, ptmLegend, psmToExport, outputFolder, 200);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR_SequenceAnnotation.png")));

            MetaDrawSettings.ExportType = "Jpeg";
            metadrawLogic.ExportSequenceCoverage(textCanvas, mapCanvas, outputFolder, psmToExport);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR_SequenceCoverage.jpeg")));
            metadrawLogic.ExportAnnotatedSequence(sequenceAnnotationCanvas, ptmLegend, psmToExport, outputFolder, 200);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR_SequenceAnnotation.jpeg")));

            MetaDrawSettings.ExportType = "Tiff";
            metadrawLogic.ExportSequenceCoverage(textCanvas, mapCanvas, outputFolder, psmToExport);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR_SequenceCoverage.tiff")));
            metadrawLogic.ExportAnnotatedSequence(sequenceAnnotationCanvas, ptmLegend, psmToExport, outputFolder, 200);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR_SequenceAnnotation.tiff")));

            MetaDrawSettings.ExportType = "Wmf";
            metadrawLogic.ExportSequenceCoverage(textCanvas, mapCanvas, outputFolder, psmToExport);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR_SequenceCoverage.wmf")));
            metadrawLogic.ExportAnnotatedSequence(sequenceAnnotationCanvas, ptmLegend, psmToExport, outputFolder, 200);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR_SequenceAnnotation.wmf")));

            MetaDrawSettings.ExportType = "Bmp";
            metadrawLogic.ExportSequenceCoverage(textCanvas, mapCanvas, outputFolder, psmToExport);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR_SequenceCoverage.bmp")));
            metadrawLogic.ExportAnnotatedSequence(sequenceAnnotationCanvas, ptmLegend, psmToExport, outputFolder, 200);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"120_QIVHDSGR_SequenceAnnotation.bmp")));

            Directory.Delete(outputFolder, true);

            // test that these methods create a directory if it is not already instantiated
            Assert.That(!Directory.Exists(outputFolder));
            metadrawLogic.ExportSequenceCoverage(textCanvas, mapCanvas, outputFolder, psmToExport);
            Assert.That(Directory.Exists(outputFolder));
            Directory.Delete(outputFolder, true);

            Assert.That(!Directory.Exists(outputFolder));
            metadrawLogic.ExportAnnotatedSequence(sequenceAnnotationCanvas, ptmLegend, psmToExport, outputFolder, 200);
            Assert.That(Directory.Exists(outputFolder));

            psm = metadrawLogic.FilteredListOfPsms[17];
            metadrawLogic.ExportSequenceCoverage(textCanvas, mapCanvas, outputFolder, psm);
            metadrawLogic.ExportAnnotatedSequence(sequenceAnnotationCanvas, ptmLegend, psm, outputFolder, 200);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"2_RGNVC[Common FixedCarbamidomet_SequenceAnnotation.bmp")));
            Assert.That(File.Exists(Path.Combine(outputFolder, @"2_RGNVC[Common FixedCarbamidomet_SequenceCoverage.bmp")));

            Directory.Delete(outputFolder, true);

        }


        [Test]
        public static void TestMetaDrawHistogramPlots()
        {
            SearchTask searchTask = new SearchTask();

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMetaDrawReadPsmFile");

            DbForTask db = new DbForTask(myDatabase, false);
            Directory.CreateDirectory(folderPath);

            searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "metadraw");
            string psmFile = Directory.GetFiles(folderPath).First(f => f.Contains("AllPSMs.psmtsv"));

            List<SpectrumMatchFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadTsv(psmFile, out var warnings);
            ObservableCollection<SpectrumMatchFromTsv> psms = new(parsedPsms);

            var psmDict = parsedPsms.GroupBy(p => p.FileNameWithoutExtension)
                .ToDictionary(p => p.Key, p => new ObservableCollection<SpectrumMatchFromTsv>(p));

            // check that fragment mass error was read in correctly
            Assert.That(Math.Round(-0.27631606125063707, 5), Is.EqualTo(Math.Round(psms[1].MatchedIons[1].MassErrorPpm, 5)));

            // check aspects of each histogram type:
            var plot = new PlotModelStat("Histogram of Precursor Masses", psms, psmDict);
            // Ensure axes are labeled correctly, and intervals are correct
            Assert.That(plot.Model.Axes.Count, Is.EqualTo(2));
            Assert.That(plot.Model.Axes[1].Title, Is.EqualTo("Count"));
            Assert.That(plot.Model.Axes[1].AbsoluteMinimum, Is.EqualTo(0));
            Assert.That(plot.Model.Axes[0].IntervalLength, Is.EqualTo(60));

            var plot2 = new PlotModelStat("Histogram of Precursor Charges", psms, psmDict);
            var series2 = plot2.Model.Series.ToList()[0];
            var items2 = (List<OxyPlot.Series.ColumnItem>)series2.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series2);
            Assert.That(items2[0].Value, Is.EqualTo(9));
            Assert.That(items2[1].Value, Is.EqualTo(1));

            var plot3 = new PlotModelStat("Histogram of Precursor PPM Errors (around 0 Da mass-difference notch only)",
                psms, psmDict);
            var series3 = plot3.Model.Series.ToList()[0];
            var items3 = (List<OxyPlot.Series.ColumnItem>)series3.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series3);
            Assert.That(items3[7].Value, Is.EqualTo(2));

            var plot4 = new PlotModelStat("Histogram of Fragment Charges",
                psms, psmDict);
            var series4 = plot4.Model.Series.ToList()[0];
            var items4 = (List<OxyPlot.Series.ColumnItem>)series4.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series4);
            Assert.That(items4[0].Value, Is.EqualTo(101));

            var plot5 = new PlotModelStat("Histogram of Precursor m/z",
                psms, psmDict);
            var series5 = plot5.Model.Series.ToList()[0];
            var items5 = (List<OxyPlot.Series.ColumnItem>)series5.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series5);
            Assert.That(items5.Count, Is.EqualTo(5));
            Assert.That(items5[0].Value, Is.EqualTo(5));

            var plot6 = new PlotModelStat("Histogram of PTM Spectral Counts",
                psms, psmDict);
            var series6 = plot6.Model.Series.ToList()[0];
            var items6 = (List<OxyPlot.Series.ColumnItem>)series6.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series6);
            Assert.That(items6.Count, Is.EqualTo(1));
            Assert.That(items6[0].Value, Is.EqualTo(2));

            var plot7 = new PlotModelStat("Precursor PPM Error vs. RT",
                psms, psmDict);
            var series7 = plot7.Model.Series.ToList()[0];
            var points7 = (List<OxyPlot.Series.ScatterPoint>)series7.GetType()
                .GetProperty("Points", BindingFlags.Public | BindingFlags.Instance).GetValue(series7);
            Assert.That(points7.Count, Is.EqualTo(9));
            Assert.That(points7[1].X, Is.EqualTo(42.07841));
            Assert.That(points7[1].Y, Is.EqualTo(-1.48));
            Assert.That(points7[1].Tag, Is.EqualTo("LSRIDTPK"));

            var plot8 = new PlotModelStat("Predicted RT vs. Observed RT",
                psms, psmDict);
            var series8 = plot8.Model.Series.ToList()[0];
            var points8 = (List<OxyPlot.Series.ScatterPoint>)series8.GetType()
                .GetProperty("Points", BindingFlags.Public | BindingFlags.Instance).GetValue(series8);
            Assert.That(points8.Count, Is.EqualTo(10));
            Assert.That(points8[7].X, Is.EqualTo(42.06171));
            Assert.That(points8[7].Y, Is.EqualTo(19.00616880619646));
            Assert.That(points8[7].Tag, Is.EqualTo("AFISYHDEAQK"));

            var plot9 = new PlotModelStat("Histogram of Fragment PPM Errors",
                psms, psmDict);
            var series9 = plot9.Model.Series.ToList()[0];
            var items9 = (List<OxyPlot.Series.ColumnItem>)series9.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series9);
            Assert.That(items9[11].Value, Is.EqualTo(18));

            var plot10 = new PlotModelStat("Histogram of Hydrophobicity scores",
                psms, psmDict);
            var series10 = plot10.Model.Series.ToList()[0];
            var items10 = (List<OxyPlot.Series.ColumnItem>)series10.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series10);
            Assert.That(items10.Count, Is.EqualTo(5));
            Assert.That(items10[1].Value, Is.EqualTo(3));

            //test variant plotting
            string variantFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\VariantCrossTest.psmtsv");
            List<string> warningsVariants = new List<string>();
            List<PsmFromTsv> parsedPsmsWithVariants;
            parsedPsmsWithVariants = SpectrumMatchTsvReader.ReadPsmTsv(variantFile, out warningsVariants);
            ObservableCollection<SpectrumMatchFromTsv> psmsWithVariants = new(parsedPsmsWithVariants);

            var psmVariantDict = psmsWithVariants.GroupBy(p => p.FileNameWithoutExtension)
                .ToDictionary(p => p.Key, p => new ObservableCollection<SpectrumMatchFromTsv>(p));

            var variantPlot1 = new PlotModelStat("Precursor PPM Error vs. RT", psmsWithVariants, psmVariantDict);
            var variantSeries1 = variantPlot1.Model.Series.ToList()[0];
            var variantPoints1 = (List<OxyPlot.Series.ScatterPoint>)variantSeries1.GetType()
                .GetProperty("Points", BindingFlags.Public | BindingFlags.Instance).GetValue(variantSeries1);
            Assert.That(variantPoints1.Count, Is.EqualTo(1));
            Assert.That(variantPoints1[0].X, Is.EqualTo(97.8357));
            Assert.That(variantPoints1[0].Y, Is.EqualTo(0.35));
            Assert.That(variantPoints1[0].Tag, Is.EqualTo("MQVDQEEPHVEEQQQQTPAENKAESEEMETSQAGSK"));

            var variantPlot2 = new PlotModelStat("Predicted RT vs. Observed RT", psmsWithVariants, psmVariantDict);
            var variantSeries2 = variantPlot2.Model.Series.ToList()[0];
            var variantPoints2 = (List<OxyPlot.Series.ScatterPoint>)variantSeries2.GetType()
                .GetProperty("Points", BindingFlags.Public | BindingFlags.Instance).GetValue(variantSeries2);
            Assert.That(variantPoints2.Count, Is.EqualTo(1));
            Assert.That(variantPoints2[0].X, Is.EqualTo(97.8357));
            Assert.That(variantPoints2[0].Y, Is.EqualTo(16.363848874371111));
            Assert.That(variantPoints2[0].Tag, Is.EqualTo("MQVDQEEPHVEEQQQQTPAENKAESEEMETSQAGSK"));

            // Histogram of Fragment Ion Types by Count
            var plot11 = new PlotModelStat("Histogram of Fragment Ion Types by Count", psms, psmDict);
            var series11 = plot11.Model.Series.ToList()[0];
            var items11 = (List<OxyPlot.Series.ColumnItem>)series11.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series11);
            Assert.That(items11.Count, Is.GreaterThan(0)); // At least one fragment type
            Assert.That(items11.Sum(i => i.Value), Is.GreaterThan(0)); // At least one count

            // Histogram of Fragment Ion Types by Intensity
            var plot12 = new PlotModelStat("Histogram of Fragment Ion Types by Intensity", psms, psmDict);
            var series12 = plot12.Model.Series.ToList()[0];
            var items12 = (List<OxyPlot.Series.ColumnItem>)series12.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series12);
            Assert.That(items12.Count, Is.GreaterThan(0)); // At least one fragment type
            Assert.That(items12.Sum(i => i.Value), Is.GreaterThan(0)); // At least one intensity

            // Histogram of Ids by Retention Time
            var plot13 = new PlotModelStat("Histogram of Ids by Retention Time", psms, psmDict);
            var series13 = plot13.Model.Series.ToList()[0];
            var items13 = (List<OxyPlot.Series.ColumnItem>)series13.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series13);
            Assert.That(items13.Count, Is.GreaterThan(0)); // At least one retention time bin
            Assert.That(items13.Sum(i => i.Value), Is.GreaterThan(0)); // At least one ID

            // Histogram of Missed Cleavages
            var plot14 = new PlotModelStat("Histogram of Missed Cleavages", psms, psmDict);
            var series14 = plot14.Model.Series.ToList()[0];
            var items14 = (List<OxyPlot.Series.ColumnItem>)series14.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series14);
            Assert.That(items14.Count, Is.GreaterThan(0)); // At least one missed cleavage bin
            Assert.That(items14.Sum(i => i.Value), Is.GreaterThanOrEqualTo(0)); // Non-negative count

            // Test normalization for missed cleavages
            MetaDrawSettings.NormalizeHistogramToFile = false;
            var plotMissedNormOff = new PlotModelStat("Histogram of Missed Cleavages", psms, psmDict);
            var seriesMissedNormOff = plotMissedNormOff.Model.Series.ToList()[0];
            var itemsMissedNormOff = (List<OxyPlot.Series.ColumnItem>)seriesMissedNormOff.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(seriesMissedNormOff);
            double sumMissedNormOff = itemsMissedNormOff.Sum(i => i.Value);

            MetaDrawSettings.NormalizeHistogramToFile = true;
            var plotMissedNormOn = new PlotModelStat("Histogram of Missed Cleavages", psms, psmDict);
            var seriesMissedNormOn = plotMissedNormOn.Model.Series.ToList()[0];
            var itemsMissedNormOn = (List<OxyPlot.Series.ColumnItem>)seriesMissedNormOn.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(seriesMissedNormOn);
            double sumMissedNormOn = itemsMissedNormOn.Sum(i => i.Value);

            Assert.That(sumMissedNormOn, Is.LessThanOrEqualTo(1.0));
            Assert.That(sumMissedNormOff, Is.GreaterThanOrEqualTo(sumMissedNormOn));

            // Reset normalization for other tests
            MetaDrawSettings.NormalizeHistogramToFile = false;

            // Test normalization OFF
            MetaDrawSettings.NormalizeHistogramToFile = false;
            var plotNormOff = new PlotModelStat("Histogram of Fragment Ion Types by Count", psms, psmDict);
            var seriesNormOff = plotNormOff.Model.Series.ToList()[0];
            var itemsNormOff = (List<OxyPlot.Series.ColumnItem>)seriesNormOff.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(seriesNormOff);
            double sumNormOff = itemsNormOff.Sum(i => i.Value);

            // Test normalization ON
            MetaDrawSettings.NormalizeHistogramToFile = true;
            var plotNormOn = new PlotModelStat("Histogram of Fragment Ion Types by Count", psms, psmDict);
            var seriesNormOn = plotNormOn.Model.Series.ToList()[0];
            var itemsNormOn = (List<OxyPlot.Series.ColumnItem>)seriesNormOn.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(seriesNormOn);
            double sumNormOn = itemsNormOn.Sum(i => i.Value);

            // When normalized, the sum should be <= 1 (since values are divided by total count)
            Assert.That(sumNormOn, Is.LessThanOrEqualTo(1.0));
            Assert.That(sumNormOff, Is.GreaterThan(sumNormOn));

            Directory.Delete(folderPath, true);
        }

        [Test]
        public static void TestCrosslinkSpectralLibraryReading()
        {
            string xlTestDataFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData");
            var interLinkResults = File.ReadAllLines(Path.Combine(xlTestDataFolder, @"XL_Interlinks.tsv"));
        }

        [Test]
        public void ExportPlot_RefragmentationWithAdditionalFragmentIons_WritesExpectedIons()
        {
            // Arrange
            var logic = new MetaDrawLogic();
            MetaDrawSettings.ExportType = "Png";
            string dataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchSingleSpectra.mzML");
            string psmFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TDGPTMDSearchResults.psmtsv");

            // Load
            logic.SpectraFilePaths.Add(dataFilePath);
            logic.SpectralMatchResultFilePaths.Add(psmFilePath);
            logic.LoadFiles(true, true);
            var psm = logic.FilteredListOfPsms.First();

            // Add to spectrumMatches
            var spectrumMatches = new List<SpectrumMatchFromTsv> { psm };

            // Set up a FragmentationReanalysisViewModel with additional ions (e.g., c, zDot)
            var reFragment = new FragmentationReanalysisViewModel();
            foreach (var frag in reFragment.PossibleProducts)
            {
                // Enable b, y, c, and zDot ions
                frag.Use = frag.ProductType == ProductType.b || frag.ProductType == ProductType.y ||
                           frag.ProductType == ProductType.bWaterLoss || frag.ProductType == ProductType.yAmmoniaLoss ||
                           frag.ProductType == ProductType.c || frag.ProductType == ProductType.zDot;
            }
            reFragment.Persist = true;

            // Set up dummy plotView and canvas
            var plotView = new PlotView { Name = "plotView" };
            var stationaryCanvas = new Canvas();
            var parentChildScanPlotsView = new ParentChildScanPlotsView();
            string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "RefragmentTest");
            Directory.CreateDirectory(tempDir);

            // Act
            // Export without refragmentation 
            string fileNoRefragment = Path.Combine(tempDir, "no_refragment.png");
            logic.ExportPlot(plotView, stationaryCanvas, spectrumMatches, parentChildScanPlotsView, tempDir, out var errorsNoRefragment, null, new System.Windows.Vector(), null);
            string exportedFileNoRefragment = Directory.GetFiles(tempDir, "*.png").FirstOrDefault();
            Assert.That(exportedFileNoRefragment, Is.Not.Null, "Exported file without refragmentation should exist.");

            // Rename File so it is not overriden by next export. 
            File.Move(exportedFileNoRefragment, fileNoRefragment);

            // Export with refragmentation
            logic.ExportPlot(plotView, stationaryCanvas, spectrumMatches, parentChildScanPlotsView, tempDir, out var errorsRefragment, null, new System.Windows.Vector(), reFragment);
            string exportedFileRefragment = Directory.GetFiles(tempDir, "*.png").FirstOrDefault();
            Assert.That(exportedFileRefragment, Is.Not.Null, "Exported file with refragmentation should exist.");

            // Assert: file with refragmentation should be larger
            var sizeNoRefragment = new FileInfo(fileNoRefragment).Length;
            var sizeRefragment = new FileInfo(exportedFileRefragment).Length;
            Assert.That(sizeRefragment, Is.GreaterThan(sizeNoRefragment), "Refragmented export should be larger due to more annotated ions.");

            // Assert: colors for c and zDot ions are present in the PNG with refragmentation but not in the other
            var cColor = MetaDrawSettings.ProductTypeToColor[ProductType.c];
            var zDotColor = MetaDrawSettings.ProductTypeToColor[ProductType.zDot];

            // Use using statements to ensure Bitmaps are disposed immediately after use
            bool noRefragmentHasC, noRefragmentHasZDot, refragmentHasC, refragmentHasZDot;
            using (var bmpNoRefragment = new Bitmap(fileNoRefragment))
            using (var bmpRefragment = new Bitmap(exportedFileRefragment))
            {
                bool HasColor(Bitmap bmp, OxyColor color)
                {
                    for (int y = 0; y < bmp.Height; y++)
                    {
                        for (int x = 0; x < bmp.Width; x++)
                        {
                            var px = bmp.GetPixel(x, y);
                            if (px.R == color.R && px.G == color.G && px.B == color.B)
                                return true;
                        }
                    }
                    return false;
                }

                noRefragmentHasC = HasColor(bmpNoRefragment, cColor);
                noRefragmentHasZDot = HasColor(bmpNoRefragment, zDotColor);
                refragmentHasC = HasColor(bmpRefragment, cColor);
                refragmentHasZDot = HasColor(bmpRefragment, zDotColor);
            }

            Assert.That(noRefragmentHasC || noRefragmentHasZDot, Is.False,
                "No c or zDot ion colors should be present in the PNG without refragmentation.");
            Assert.That(refragmentHasC || refragmentHasZDot, Is.True,
                "c or zDot ion colors should be present in the PNG with refragmentation.");

            Assert.That(errorsNoRefragment, Is.Null, "No errors should be reported for no refragmentation.");
            Assert.That(errorsRefragment, Is.Null, "No errors should be reported for refragmentation.");

            // Clean up
            Directory.Delete(tempDir, true);
        }

        [Test] // Ensures no plot names accidently get deleted. 
        public static void PlotNamesDoNotChange()
        {
            var plotNames = PlotModelStat.PlotNames;
            Assert.That(plotNames.Count, Is.EqualTo(14));
        }
    }
}