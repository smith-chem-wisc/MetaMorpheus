using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using TaskLayer;
using EngineLayer;
using System.Linq;
using EngineLayer.DatabaseLoading;
using NUnit.Framework.Legacy;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;

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
            Assert.That(sections.Length, Is.EqualTo(48));

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

            if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
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
            Assert.That(sections.Length, Is.EqualTo(48));

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

            Directory.Delete(outputFolder, true);
        }
    }
}
