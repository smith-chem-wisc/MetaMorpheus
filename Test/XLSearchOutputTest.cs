using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using TaskLayer;
using EngineLayer;

namespace Test
{
    [TestFixture]
    public static class XLSearchOutputTest
    {
        [Test]
        public static void WriteTsvTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlOutputTest1");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA_DSS_23747.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA.fasta");

            Directory.CreateDirectory(outputFolder);

            XLSearchTask xLSearch = new XLSearchTask();
            xLSearch.XlSearchParameters.CrosslinkAtCleavageSite = true;
            xLSearch.XlSearchParameters.Crosslinker = new Crosslinker("K", "K", "DSS", false, "HCD", 138.06808, 0, 0, 138.06808, 156.0786, 155.0946, 259.142);
            xLSearch.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");

            var resultsPath = File.ReadAllLines(Path.Combine(outputFolder, @"XL_Intralinks.tsv"));
            var sections = resultsPath[1].Split('\t');
            Assert.That(resultsPath.Length == 3);
            Assert.That(sections.Length == 43);

            var resultsPath_Deadend = File.ReadAllLines(Path.Combine(outputFolder, @"Deadends.tsv"));
            Assert.That(resultsPath_Deadend.Length == 2);

            Directory.Delete(outputFolder, true);
        }
    }
}
