using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.Collections.Generic;
using System.IO;
using TaskLayer;
using EngineLayer;
using System.Linq;
using NUnit.Framework.Legacy;
using Omics.Fragmentation;

namespace Test
{
    [TestFixture]
    public static class XLSearchOutputTest
    {
        [Test]
        public static void WriteTsvTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlOutputTestFile");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\2017-11-21_XL_DSSO_Ribosome_RT60min_28800-28898.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\RibosomeGO.fasta");

            Directory.CreateDirectory(outputFolder);

            XLSearchTask xLSearch = new XLSearchTask();
            xLSearch.XlSearchParameters.CrosslinkAtCleavageSite = true;
            xLSearch.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");

            var resultsPath = File.ReadAllLines(Path.Combine(outputFolder, @"XL_Intralinks.tsv"));
            var sections = resultsPath[1].Split('\t');
            Assert.That(resultsPath.Length > 1);
            Assert.AreEqual(sections.Length, 48);

            var resultsPath_Inter = File.ReadAllLines(Path.Combine(outputFolder, @"XL_Interlinks.tsv"));
            Assert.That(resultsPath_Inter.Length > 1);

            var resultsPath_Deadend = File.ReadAllLines(Path.Combine(outputFolder, @"Deadends.tsv"));
            Assert.That(resultsPath_Deadend.Length >1);

            var resultsPath_loop = File.ReadAllLines(Path.Combine(outputFolder, @"Looplinks.tsv"));
            Assert.That(resultsPath_loop.Length >1);

            var resultsPath_single = File.ReadAllLines(Path.Combine(outputFolder, @"SinglePeptides.tsv"));
            Assert.That(resultsPath_single.Length >1);

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void WriteXlSpectralLibraryTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlOutputTestFile");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\2017-11-21_XL_DSSO_Ribosome_RT60min_28800-28898.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\RibosomeGO.fasta");

            Directory.CreateDirectory(outputFolder);

            XLSearchTask xLSearch = new XLSearchTask
            {
                XlSearchParameters = new XlSearchParameters
                {
                    WriteSpectralLibrary = true,
                    CrosslinkAtCleavageSite = true
                }
            };
            xLSearch.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");

            var resultsPath = File.ReadAllLines(Path.Combine(outputFolder, @"XL_Intralinks.tsv"));
            var sections = resultsPath[1].Split('\t');
            Assert.That(resultsPath.Length > 1);
            Assert.AreEqual(sections.Length, 48);

            var resultsPath_Inter = File.ReadAllLines(Path.Combine(outputFolder, @"XL_Interlinks.tsv"));
            Assert.That(resultsPath_Inter.Length > 1);

            var resultsPath_Deadend = File.ReadAllLines(Path.Combine(outputFolder, @"Deadends.tsv"));
            Assert.That(resultsPath_Deadend.Length > 1);

            var resultsPath_loop = File.ReadAllLines(Path.Combine(outputFolder, @"Looplinks.tsv"));
            Assert.That(resultsPath_loop.Length > 1);

            var resultsPath_single = File.ReadAllLines(Path.Combine(outputFolder, @"SinglePeptides.tsv"));
            Assert.That(resultsPath_single.Length > 1);

            var list = Directory.GetFiles(outputFolder, "*.*", SearchOption.AllDirectories);
            string matchingvalue = list.Where(p => p.Contains("SpectralLibrary")).First().ToString();
            var testLibraryWithoutDecoy = new SpectralLibrary(new List<string> { Path.Combine(outputFolder, matchingvalue) });

            // Check that SinglePeptides are written
            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("GNVINLSLGFSHPVDHQLPAGITAEC[Common Fixed:Carbamidomethyl on C]PTQTEIVLK", 4, out var spectrum));
            Product productWithNeutralLoss =
                new Product(ProductType.Y, FragmentationTerminus.C, 100, 1, 1, neutralLoss: 10.0);
            Product productWithNeutralLoss20 =
                new Product(ProductType.Y, FragmentationTerminus.C, 100, 1, 1, neutralLoss: 20.0);
            // Check that interLinks
            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("GVTVDKMTELR(6)SFTFVTKTPPAAVLLK(7)", 2, out var interLinkSpectrum));
            Assert.That(interLinkSpectrum.MatchedFragmentIons.Count, Is.EqualTo(17));
            Assert.That(interLinkSpectrum is CrosslinkLibrarySpectrum);
            CrosslinkLibrarySpectrum interSpectrum = (CrosslinkLibrarySpectrum)interLinkSpectrum;
            Assert.That(interSpectrum.BetaPeptideSpectrum.MatchedFragmentIons.Count, Is.EqualTo(13));
            Assert.AreEqual(interSpectrum.AlphaPeptideSequence, "GVTVDKMTELR");
            Assert.AreEqual(interSpectrum.BetaPeptideSequence, "SFTFVTKTPPAAVLLK");
            Assert.True(interSpectrum.BetaPeptideSpectrum.IsBetaPeptide);
            interLinkSpectrum.MatchedFragmentIons.Add(new MatchedFragmentIon(productWithNeutralLoss, 100, 100, 1));
            CrosslinkLibrarySpectrum spectrumDup = (CrosslinkLibrarySpectrum)interLinkSpectrum;
            spectrumDup.BetaPeptideSpectrum.MatchedFragmentIons.Add(new MatchedFragmentIon(productWithNeutralLoss20, 100, 100, 1));
            var spectrumString = spectrumDup.ToString();
            // Check neutral loss fragments are written correctly
            StringAssert.Contains("\"Y1^1-10/0ppm\"", spectrumString);
            StringAssert.Contains("\"Y1^1-20/0ppm\"\tBetaPeptideIon", spectrumString);
            

            // Check intraLinks
            Assert.That(testLibraryWithoutDecoy.TryGetSpectrum("LLDNAAADLAAISGQKPLITKAR(21)ITLNMGVGEAIADKK(14)", 5, out var intraLinkSpectrum));
            
            Assert.That(intraLinkSpectrum is CrosslinkLibrarySpectrum);
            CrosslinkLibrarySpectrum intraSpectrum = (CrosslinkLibrarySpectrum)interLinkSpectrum;
            Assert.That(intraSpectrum.BetaPeptideSpectrum.MatchedFragmentIons.Count, Is.GreaterThan(0));

            testLibraryWithoutDecoy.CloseConnections();

            Directory.Delete(outputFolder, true);
        }
    }
}
