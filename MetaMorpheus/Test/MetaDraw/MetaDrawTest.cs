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
using MassSpectrometry;
using NUnit.Framework;
using OxyPlot.Series;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Readers;
using TaskLayer;

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

            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);

            Assert.AreEqual(10, parsedPsms.Count);
            Assert.AreEqual(0, warnings.Count);

            Directory.Delete(folderPath, true);
        }


        [Test]
        public static void TestMetaDrawReadPsmFileWithoutSpectralAngle()
        {
            //test if the reader still works when psm file doesn't contain spectral angle as header.
            string noSA = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\noSAreaderTest.psmtsv");
            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(noSA, out var warnings);


            Assert.AreEqual(15, parsedPsms.Count);
            Assert.AreEqual(0, warnings.Count);
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

            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);

            Assert.AreEqual(-1, parsedPsms.First().SpectralAngle);
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
        public static void MetaDraw_TestStationarySequencePositioning()
        {
            var metadrawLogic = new MetaDrawLogic();
            metadrawLogic.PsmResultFilePaths.Add(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TopDownTestData\TDGPTMDSearchResults.psmtsv"));
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

                List<MatchedFragmentIon> matchedIons = psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.AminoAcidPosition > MetaDrawSettings.FirstAAonScreenIndex &&
                                                       p.NeutralTheoreticalProduct.AminoAcidPosition < (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen)).ToList();
                int psmStartResidue = int.Parse(psm.StartAndEndResiduesInProtein.Split("to")[0].Replace("[", ""));
                var startAA = (MetaDrawSettings.FirstAAonScreenIndex + psmStartResidue).ToString();
                var endAA = (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen + psmStartResidue - 1).ToString();


                var expected = modifiedBaseSeq.Length + matchedIons.Count + fullSequence.Count(p => p == '[') + startAA.Length + endAA.Length + 2;
                Assert.AreEqual(metadrawLogic.StationarySequence.SequenceDrawingCanvas.Children.Count, expected);
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
            List<MatchedFragmentIon> matchedIons = psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.AminoAcidPosition > MetaDrawSettings.FirstAAonScreenIndex &&
                                                   p.NeutralTheoreticalProduct.AminoAcidPosition < (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen)).ToList();
            Assert.That(metadrawLogic.StationarySequence.SequenceDrawingCanvas.Children.Count == modifiedBaseSeq.Length + matchedIons.Count + fullSequence.Count(p => p == '[') + 
                psm.StartAndEndResiduesInProtein.Replace("[","").Replace("]","").Replace("to", "").Replace(" ", "").Length + 2);

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
            matchedIons = modPsm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.AminoAcidPosition > MetaDrawSettings.FirstAAonScreenIndex &&
                                                p.NeutralTheoreticalProduct.AminoAcidPosition < (MetaDrawSettings.FirstAAonScreenIndex + MetaDrawSettings.NumberOfAAOnScreen)).ToList();
            Assert.That(metadrawLogic.StationarySequence.SequenceDrawingCanvas.Children.Count == modifiedBaseSeq.Length + matchedIons.Count + fullSequence.Count(p => p == '[') + 
                psm.StartAndEndResiduesInProtein.Replace("[", "").Replace("]", "").Replace("to", "").Replace(" ", "").Length + 2);


            // get scan from psm
            var scan = metadrawLogic.GetMs2ScanFromPsm(psm);
            Assert.That(scan.OneBasedScanNumber, Is.EqualTo(psm.Ms2ScanNumber));
            Assert.That(scan.MsnOrder, Is.EqualTo(2));

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
            Assert.That(!metadrawLogic.PsmResultFilePaths.Any());
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
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
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
            Assert.AreEqual(mirrorPlotPeaks.Count, 59);

            var plotAxes = plotView.Model.Axes;
            Assert.That(plotAxes.Count == 2);

            // write pdf
            var psmsToExport = metadrawLogic.FilteredListOfPsms.Where(p => p.UniqueSequence == "LLDNAAADLAAISGQKPLITKAR(21)ITLNMGVGEAIADKK(14)").Take(1).ToList();
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
            Assert.That(!metadrawLogic.PsmResultFilePaths.Any());
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
            Assert.Greater(c, 0);

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
            Assert.Greater(c, 0);

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
            Assert.AreEqual(2, plotAxes.Count);

            // test that base sequence annotation was drawn
            Assert.Greater(metadrawLogic.StationarySequence.SequenceDrawingCanvas.Children.Count, 0);

            // test that the plots were drawn in the parent/child view
            Assert.AreEqual(2, parentChildView.Plots.Count);

            // test parent scan
            var parentPlot = parentChildView.Plots[0];
            Assert.AreEqual("Scan: 27 Dissociation Type: HCD MsOrder: 2 Selected Mz: 924.12 Retention Time: 32.65", parentPlot.SpectrumLabel);
            int numAnnotatedResidues = psm.BaseSeq.Length;
            int numAnnotatedIons = psm.MatchedIons.Count(p => p.NeutralTheoreticalProduct.ProductType != ProductType.M
                && p.NeutralTheoreticalProduct.ProductType != ProductType.D);
            int numAnnotatedMods = psm.FullSequence.Count(p => p == '[');
            Assert.AreEqual(numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods, parentPlot.TheCanvas.Children.Count);

            peak = (LineSeries)parentPlot.Plot.Model.Series[0]; // the first m/z peak
            peakPoints = peak.Points;
            Assert.AreEqual(101.07, Math.Round(peakPoints[0].X, 2)); // m/z
            Assert.AreEqual(101.07, (Math.Round(peakPoints[1].X, 2)));
            Assert.AreEqual(0, (int)peakPoints[0].Y); // intensity
            Assert.AreEqual(3847, (int)peakPoints[1].Y);

            // test child scan
            var childPlot = parentChildView.Plots[1];
            Assert.AreEqual("Scan: 30 Dissociation Type: EThcD MsOrder: 2 Selected Mz: 924.12 RetentionTime: 32.66",
                childPlot.SpectrumLabel);
            Assert.Greater(childPlot.TheCanvas.Children.Count, 0);
            numAnnotatedResidues = psm.BaseSeq.Length;
            numAnnotatedIons = psm.ChildScanMatchedIons[30]
                .Count(p => p.NeutralTheoreticalProduct.ProductType != ProductType.M
                && p.NeutralTheoreticalProduct.ProductType != ProductType.D);
            numAnnotatedMods = psm.FullSequence.Count(p => p == '[');
            Assert.AreEqual(numAnnotatedResidues + numAnnotatedIons + numAnnotatedMods, childPlot.TheCanvas.Children.Count);

            peak = (LineSeries)childPlot.Plot.Model.Series[0]; // the first m/z peak
            peakPoints = peak.Points;
            Assert.AreEqual(126.06, Math.Round(peakPoints[0].X, 2)); // m/z
            Assert.AreEqual(126.06, Math.Round(peakPoints[1].X, 2));
            Assert.AreEqual(0, (int)peakPoints[0].Y); // intensity
            Assert.AreEqual(8496, (int)peakPoints[1].Y);

            // write pdf
            var psmsToExport = metadrawLogic.FilteredListOfPsms.Where(p => p.FullSequence == "STTAVQTPTSGEPLVST[O-Glycosylation:H1N1 on X]SEPLSSK").ToList();
            metadrawLogic.ExportPlot(plotView, metadrawLogic.StationarySequence.SequenceDrawingCanvas, psmsToExport, parentChildView, outputFolder, out errors);

            // test that pdf exists
            Assert.That(File.Exists(Path.Combine(outputFolder, @"27_STTAVQTPTSGEPLVST[O-Glycosylat.pdf"))); // parent scan
            Assert.That(File.Exists(Path.Combine(outputFolder, @"30_STTAVQTPTSGEPLVST[O-Glycosylat.pdf"))); // child scan
            Directory.Delete(outputFolder, true);

            Canvas ptmLegend = new();
            System.Windows.Size legendSize = new(100, 100);
            ptmLegend.Measure(legendSize);
            ptmLegend.Arrange(new Rect(legendSize));
            ptmLegend.UpdateLayout();   
            Vector ptmLegendVector = new(10, 10);
            metadrawLogic.ExportPlot(plotView, metadrawLogic.StationarySequence.SequenceDrawingCanvas, psmsToExport, parentChildView, outputFolder, out errors, ptmLegend, ptmLegendVector);
            Assert.That(File.Exists(Path.Combine(outputFolder, @"27_STTAVQTPTSGEPLVST[O-Glycosylat.pdf"))); // parent scan
            Assert.That(File.Exists(Path.Combine(outputFolder, @"30_STTAVQTPTSGEPLVST[O-Glycosylat.pdf"))); // child scan

            // test that a directory is created if it does not exist
            Directory.Delete(outputFolder, true);
            metadrawLogic.ExportPlot(plotView, metadrawLogic.StationarySequence.SequenceDrawingCanvas, psmsToExport, parentChildView, outputFolder, out errors);
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
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(true, true);
            Assert.That(!errors.Any());
            Assert.That(metadrawLogic.FilteredListOfPsms.Any());

            // fix the scan number due to the trimmed spectra file
            foreach (var psm in metadrawLogic.FilteredListOfPsms)
            {
                var type = psm.GetType();
                var field = type.GetField("<Ms2ScanNumber>k__BackingField", BindingFlags.Instance | BindingFlags.NonPublic);
                field.SetValue(psm, psm.Ms2ScanNumber + 27300);
            }
            var metaDrawDynamicScanConnection = (Dictionary<string, MsDataFile>)metadrawLogic?.GetType()
                .GetField("MsDataFiles", BindingFlags.Instance | BindingFlags.NonPublic)
                ?.GetValue(metadrawLogic);


            metadrawLogic.FilterPsmsToChimerasOnly();
            // test plotting on each instance of chimeras in this dataset
            var plotView = new OxyPlot.Wpf.PlotView() { Name = "chimeraPlot" };
            foreach (var chimeraGroup in metadrawLogic.FilteredListOfPsms
                         .GroupBy(p => p.Ms2ScanNumber))
            {
                Assert.That(chimeraGroup.Count(), Is.GreaterThanOrEqualTo(2));
                MsDataScan chimericScan = metaDrawDynamicScanConnection.First().Value
                    .GetOneBasedScanFromDynamicConnection(chimeraGroup.First().Ms2ScanNumber);

                // plot the first chimera and test the results
                metadrawLogic.DisplayChimeraSpectra(plotView, chimeraGroup.ToList(), out errors);
                Assert.That(errors == null || !errors.Any());

                // test plot was drawn
                var model = plotView.Model;
                Assert.That(model, Is.Not.Null);
                Assert.That(model.Equals(metadrawLogic.ChimeraSpectrumMatchPlot.Model));
                Assert.That(plotView.Model.Axes.Count == 2);

                var peakPoints = ((LineSeries)model.Series[0]).Points;
                
                Assert.That(Math.Round(peakPoints[0].X, 2) == Math.Round(chimericScan.MassSpectrum.XArray[0], 2)); // m/z
                Assert.That(Math.Round(peakPoints[1].X, 2) == Math.Round(chimericScan.MassSpectrum.XArray[0], 2));
                Assert.That((int)peakPoints[0].Y == 0); // intensity
                Assert.That((int)peakPoints[1].Y == (int)chimericScan.MassSpectrum.YArray[0]);

                // all matched ions were drawn
                int drawnIonsNotDefaultColor = model.Series.Count(p => ((LineSeries)p).Color != MetaDrawSettings.UnannotatedPeakColor);
                List<MatchedFragmentIon> fragments = new();
                chimeraGroup.Select(p => p.MatchedIons).ForEach(m => fragments.AddRange(m));
                Assert.That(drawnIonsNotDefaultColor, Is.EqualTo(fragments.Count));

                // shared matched ions are default color
                int drawnIonsShared = model.Series.Count(p => ((LineSeries)p).Color == ChimeraSpectrumMatchPlot.MultipleProteinSharedColor);
                var sharedIons = fragments.GroupBy(p => p)
                    .Where(m => m.Count() > 1)
                    .SelectMany(n => n).ToList();
                if (sharedIons.Any() || drawnIonsShared >= 1)
                {
                    int distinctMatchedSharedIons = sharedIons.Distinct().Count();
                    Assert.That(sharedIons.Count - distinctMatchedSharedIons == drawnIonsShared);
                }

                // unshared peaks are the correct color
                var unsharedIons = fragments.GroupBy(p => p)
                    .Where(m => m.Count() == 1)
                    .SelectMany(n => n).ToList();
                for (var i = 0; i < chimeraGroup.Count(); i++)
                {
                    var chimera = chimeraGroup.ElementAt(i);
                    var chimeraSpecificPeaks = unsharedIons.Intersect(chimera.MatchedIons).ToList();
                    var chimeraSharedPeaks = sharedIons.Intersect(chimera.MatchedIons).ToList();
                    int drawnIonsOfSpecificID = model.Series.Count(p => ChimeraSpectrumMatchPlot.ColorByProteinDictionary[i].Any(m => m == ((LineSeries)p).Color));

                    if (i == 0)
                        Assert.That(drawnIonsOfSpecificID - chimeraSharedPeaks.Count == chimeraSpecificPeaks.Count);
                    else
                        Assert.That(drawnIonsOfSpecificID == chimeraSpecificPeaks.Count);
                }

                // test with different drawing settings
                MetaDrawSettings.AnnotateCharges = true;
                MetaDrawSettings.AnnotateMzValues = true;
                metadrawLogic.DisplayChimeraSpectra(plotView, chimeraGroup.ToList(), out errors);
                Assert.That(errors == null || !errors.Any());

                MetaDrawSettings.DisplayIonAnnotations = false;
                metadrawLogic.DisplayChimeraSpectra(plotView, chimeraGroup.ToList(), out errors);
                Assert.That(errors == null || !errors.Any());

                MetaDrawSettings.DisplayInternalIons = false;
                metadrawLogic.DisplayChimeraSpectra(plotView, chimeraGroup.ToList(), out errors);
                Assert.That(errors == null || !errors.Any());
                MetaDrawSettings.DisplayInternalIons = true;
            }

            // test export of singlular plot
            List<PsmFromTsv> firstChimeraGroup = metadrawLogic.FilteredListOfPsms
                .GroupBy(p => p.Ms2ScanNumber).First().ToList();
            metadrawLogic.DisplayChimeraSpectra(plotView, firstChimeraGroup, out errors);
            Assert.That(errors == null || !errors.Any());
            foreach (var exportType in MetaDrawSettings.ExportTypes)
            {
                MetaDrawSettings.ExportType = exportType;
                metadrawLogic.ExportPlot(plotView, null, new List<PsmFromTsv>() { firstChimeraGroup.First() }, null, outputFolder, out errors);
                Assert.That(errors == null || !errors.Any());
                string sequence = illegalInFileName.Replace(firstChimeraGroup.First().FullSequence, string.Empty);
                string filePathWithoutDirectory = firstChimeraGroup.First().Ms2ScanNumber + "_" 
                    + (sequence.Length > 30 ? sequence.Substring(0, 30) : sequence)
                    + "." + exportType;
                Assert.That(File.Exists(Path.Combine(outputFolder, filePathWithoutDirectory)));
            }

            string export = MetaDrawSettings.ExportType;
            Canvas ptmLegend = new();
            System.Windows.Size legendSize = new(100, 100);
            ptmLegend.Measure(legendSize);
            ptmLegend.Arrange(new Rect(legendSize));
            ptmLegend.UpdateLayout();
            Vector ptmLegendVector = new(10, 10);
            metadrawLogic.ExportPlot(plotView, null, new List<PsmFromTsv>() { firstChimeraGroup.First() }, null,
                outputFolder, out errors, ptmLegend, ptmLegendVector);
            Assert.That(errors == null || !errors.Any());
            string sequenceSeq = illegalInFileName.Replace(firstChimeraGroup.First().FullSequence, string.Empty);
            string fileNameWithoutDirectory = firstChimeraGroup.First().Ms2ScanNumber + "_"
                + (sequenceSeq.Length > 30 ? sequenceSeq.Substring(0, 30) : sequenceSeq)
                + "." + export;
            Assert.That(File.Exists(Path.Combine(outputFolder, fileNameWithoutDirectory)));

            // test export of multiple plots
            List<PsmFromTsv> secondChimeraGroup = metadrawLogic.FilteredListOfPsms
                .GroupBy(p => p.Ms2ScanNumber).ToList()[1].ToList();
            metadrawLogic.DisplayChimeraSpectra(plotView, secondChimeraGroup, out errors);
            Assert.That(errors == null || !errors.Any());
            foreach (var exportType in MetaDrawSettings.ExportTypes)
            {
                MetaDrawSettings.ExportType = exportType;
                metadrawLogic.ExportPlot(plotView, null, secondChimeraGroup, null, outputFolder, out errors);
                Assert.That(errors == null || !errors.Any());

                foreach (var chimera in secondChimeraGroup)
                {
                    string sequence = illegalInFileName.Replace(chimera.FullSequence, string.Empty);
                    string filePathWithoutDirectory = chimera.Ms2ScanNumber + "_" 
                        + (sequence.Length > 30 ? sequence.Substring(0, 30) : sequence) 
                        + "." + exportType;
                    Assert.That(File.Exists(Path.Combine(outputFolder, filePathWithoutDirectory)));
                }
            }

            // test error
            metadrawLogic.CleanUpResources();
            metadrawLogic.DisplayChimeraSpectra(plotView, secondChimeraGroup, out errors);
            Assert.That(errors != null && errors.First().Equals("The spectra file could not be found for this PSM: TaGe_SA_HeLa_04_subset_longestSeq"));
            
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
            var scrollableCanvas = new Canvas();
            var parentChildScanPlotsView = new ParentChildScanPlotsView();

            // plotting PSM should produce an error because spectra are not loaded
            metadrawLogic.DisplaySpectrumMatch(plotView, psmsFromTsv.First(), parentChildScanPlotsView, out errors);
            Assert.That(errors.Any());

            // export to PDF should produce an error because spectra are not loaded
            MetaDrawSettings.NumberOfAAOnScreen = psmsFromTsv.First().BaseSeq.Length - 1;
            metadrawLogic.ExportPlot(plotView, canvas, new List<PsmFromTsv> { psmsFromTsv.First() }, parentChildScanPlotsView, outputFolder, out errors);
            Assert.That(errors.Any());

            // clean up resources
            metadrawLogic.CleanUpResources();

            // delete output
            Directory.Delete(outputFolder, true);
        }

        [Test]
        [NonParallelizable]
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
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
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
            metadrawLogic.ExportPlot(plotView, canvas, new List<PsmFromTsv> { metadrawLogic.FilteredListOfPsms.First() }, parentChildScanPlotsView, outputFolder, out errors);
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
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
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
                Assert.That(mz.ToMass(charge) - parsedIon.MassErrorDa == parsedIon.NeutralTheoreticalProduct.NeutralMass);
                Assert.That(charge == parsedIon.Charge);
                Assert.That(ion == parsedIon.NeutralTheoreticalProduct.ProductType.ToString() + parsedIon.NeutralTheoreticalProduct.FragmentNumber);
            }

            //check that the proteoform classification was correct
            Assert.IsTrue(metadrawLogic.FilteredListOfPsms[0].AmbiguityLevel.Equals("1"));

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
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
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
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
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
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
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
            Assert.That(!metadrawLogic.PsmResultFilePaths.Any());

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
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
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
            var psmsToExport = new List<PsmFromTsv> { metadrawLogic.FilteredListOfPsms.Where(p => p.FullSequence == "QIVHDSGR").First() };
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
            var psmOfInterest = PsmTsvReader.ReadTsv(psmsPath, out var warnings).First(p => p.FullSequence.Contains("LLSTEGSQ"));
            Assert.That(!warnings.Any());

            MetaDrawLogic metaDrawLogic = new();
            metaDrawLogic.PsmResultFilePaths.Add(psmsPath);
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

            psm = metadrawLogic.FilteredListOfPsms[19];
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

            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            ObservableCollection<PsmFromTsv> psms = new(parsedPsms);

            var psmDict = parsedPsms.GroupBy(p => p.FileNameWithoutExtension)
                .ToDictionary(p => p.Key, p => new ObservableCollection<PsmFromTsv>(p));

            // check that fragment mass error was read in correctly
            Assert.AreEqual(Math.Round(-0.27631606125063707, 5), Math.Round(psms[1].MatchedIons[1].MassErrorPpm, 5));
            
            // check aspects of each histogram type:
            var plot = new PlotModelStat("Histogram of Precursor Masses", psms, psmDict);
            // Ensure axes are labeled correctly, and intervals are correct
            Assert.AreEqual(2, plot.Model.Axes.Count);
            Assert.AreEqual("Count", plot.Model.Axes[1].Title);
            Assert.AreEqual(0, plot.Model.Axes[1].AbsoluteMinimum);
            Assert.AreEqual(60, plot.Model.Axes[0].IntervalLength);

            var plot2 = new PlotModelStat("Histogram of Precursor Charges", psms, psmDict);
            var series2 = plot2.Model.Series.ToList()[0];
            var items2 = (List<OxyPlot.Series.ColumnItem>)series2.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series2);
            Assert.AreEqual(items2[0].Value, 9);
            Assert.AreEqual(items2[1].Value, 1);

            var plot3 = new PlotModelStat("Histogram of Precursor PPM Errors (around 0 Da mass-difference notch only)",
                psms, psmDict);
            var series3 = plot3.Model.Series.ToList()[0];
            var items3 = (List<OxyPlot.Series.ColumnItem>)series3.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series3);
            Assert.AreEqual(items3[7].Value, 2);

            var plot4 = new PlotModelStat("Histogram of Fragment Charges",
                psms, psmDict);
            var series4 = plot4.Model.Series.ToList()[0];
            var items4 = (List<OxyPlot.Series.ColumnItem>)series4.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series4);
            Assert.AreEqual(items4[0].Value, 101);

            var plot5 = new PlotModelStat("Histogram of Precursor m/z",
                psms, psmDict);
            var series5 = plot5.Model.Series.ToList()[0];
            var items5 = (List<OxyPlot.Series.ColumnItem>)series5.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series5);
            Assert.AreEqual(items5.Count, 5);
            Assert.AreEqual(items5[0].Value, 5);

            var plot6 = new PlotModelStat("Histogram of PTM Spectral Counts",
                psms, psmDict);
            var series6 = plot6.Model.Series.ToList()[0];
            var items6 = (List<OxyPlot.Series.ColumnItem>)series6.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series6);
            Assert.AreEqual(items6.Count, 1);
            Assert.AreEqual(items6[0].Value, 2);

            var plot7 = new PlotModelStat("Precursor PPM Error vs. RT",
                psms, psmDict);
            var series7 = plot7.Model.Series.ToList()[0];
            var points7 = (List<OxyPlot.Series.ScatterPoint>)series7.GetType()
                .GetProperty("Points", BindingFlags.Public | BindingFlags.Instance).GetValue(series7);
            Assert.AreEqual(points7.Count, 9);
            Assert.AreEqual(points7[1].X, 42.07841);
            Assert.AreEqual(points7[1].Y, -1.48);
            Assert.AreEqual(points7[1].Tag, "LSRIDTPK");

            var plot8 = new PlotModelStat("Predicted RT vs. Observed RT",
                psms, psmDict);
            var series8 = plot8.Model.Series.ToList()[0];
            var points8 = (List<OxyPlot.Series.ScatterPoint>)series8.GetType()
                .GetProperty("Points", BindingFlags.Public | BindingFlags.Instance).GetValue(series8);
            Assert.AreEqual(points8.Count, 10);
            Assert.AreEqual(points8[7].X, 42.06171);
            Assert.AreEqual(points8[7].Y, 19.00616880619646);
            Assert.AreEqual(points8[7].Tag, "AFISYHDEAQK");

            var plot9 = new PlotModelStat("Histogram of Fragment PPM Errors",
                psms, psmDict);
            var series9 = plot9.Model.Series.ToList()[0];
            var items9 = (List<OxyPlot.Series.ColumnItem>)series9.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series9);
            Assert.AreEqual(items9[11].Value, 18);

            var plot10 = new PlotModelStat("Histogram of Hydrophobicity scores",
                psms, psmDict);
            var series10 = plot10.Model.Series.ToList()[0];
            var items10 = (List<OxyPlot.Series.ColumnItem>)series10.GetType()
                .GetProperty("Items", BindingFlags.Public | BindingFlags.Instance).GetValue(series10);
            Assert.AreEqual(items10.Count, 5);
            Assert.AreEqual(items10[1].Value, 3);

            //test variant plotting
            string variantFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\VariantCrossTest.psmtsv");
            List<string> warningsVariants = new List<string>();
            List<PsmFromTsv> parsedPsmsWithVariants;
            parsedPsmsWithVariants = PsmTsvReader.ReadTsv(variantFile, out warningsVariants);
            ObservableCollection<PsmFromTsv> psmsWithVariants = new(parsedPsmsWithVariants);

            var psmVariantDict = psmsWithVariants.GroupBy(p => p.FileNameWithoutExtension)
                .ToDictionary(p => p.Key, p => new ObservableCollection<PsmFromTsv>(p));

            var variantPlot1 = new PlotModelStat("Precursor PPM Error vs. RT", psmsWithVariants, psmVariantDict);
            var variantSeries1 = variantPlot1.Model.Series.ToList()[0];
            var variantPoints1 = (List<OxyPlot.Series.ScatterPoint>)variantSeries1.GetType()
                .GetProperty("Points", BindingFlags.Public | BindingFlags.Instance).GetValue(variantSeries1);
            Assert.AreEqual(variantPoints1.Count, 1);
            Assert.AreEqual(variantPoints1[0].X, 97.8357);
            Assert.AreEqual(variantPoints1[0].Y, 0.35);
            Assert.AreEqual(variantPoints1[0].Tag, "MQVDQEEPHVEEQQQQTPAENKAESEEMETSQAGSK");

            var variantPlot2 = new PlotModelStat("Predicted RT vs. Observed RT", psmsWithVariants, psmVariantDict);
            var variantSeries2 = variantPlot2.Model.Series.ToList()[0];
            var variantPoints2 = (List<OxyPlot.Series.ScatterPoint>)variantSeries2.GetType()
                .GetProperty("Points", BindingFlags.Public | BindingFlags.Instance).GetValue(variantSeries2);
            Assert.AreEqual(variantPoints2.Count, 1);
            Assert.AreEqual(variantPoints2[0].X, 97.8357);
            Assert.AreEqual(variantPoints2[0].Y, 16.363848874371111);
            Assert.AreEqual(variantPoints2[0].Tag, "MQVDQEEPHVEEQQQQTPAENKAESEEMETSQAGSK");


            Directory.Delete(folderPath, true);
        }

        [Test]
        public static void TestCrosslinkSpectralLibraryReading()
        {
            string xlTestDataFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData");
            var interLinkResults = File.ReadAllLines(Path.Combine(xlTestDataFolder, @"XL_Interlinks.tsv"));
        }
        
    }
}