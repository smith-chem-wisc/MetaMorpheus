using EngineLayer;
using MetaMorpheusGUI;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
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
            string psmFile = Directory.GetFiles(folderPath).First(p => p.Contains("AllPSMs.psmtsv"));

            List<MetaDrawPsm> parsedPsms = TsvResultReader.ReadTsv(psmFile, out var warnings);

            Assert.AreEqual(11, parsedPsms.Count);
            Assert.AreEqual(0, warnings.Count);

            Directory.Delete(folderPath, true);
        }
    }
}
