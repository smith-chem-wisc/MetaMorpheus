using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using TaskLayer;

namespace Test
{
    [TestFixture]
    class XLSearchOutputTest
    {
        [Test]
        public void WriteTsvTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XLTestData\OutputTest1");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XLTestData\BSA_DSS_23747.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XLTestData\BSA.fasta");

            XLSearchTask xLSearch = new XLSearchTask();
            xLSearch.RunTask(outputFolder, new List<DbForTask> { new DbForTask(myDatabase, false) }, new List<string> { myFile }, "test");

            var resultsPath = File.ReadAllLines(Path.Combine(outputFolder, @"XL_Interlinks.tsv"));
            var sections = resultsPath[1].Split('\t');
            Assert.That(resultsPath.Length == 2);
            Assert.That(sections.Length == 31);

        }
    }
}
