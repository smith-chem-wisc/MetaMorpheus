using EngineLayer;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using Proteomics.Fragmentation;

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
            string psmFile = Directory.GetFiles(folderPath).First(f => f.Contains("AllPSMs.psmtsv"));

            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);

            Assert.AreEqual(11, parsedPsms.Count);
            Assert.AreEqual(0, warnings.Count);

            Directory.Delete(folderPath, true);
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
        public static void TestLabelingVariantCrossingIons()
        {
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\VariantCrossTest.psmtsv");
            List<string> warnings = new List<string>();
            List<PsmFromTsv> psms;

            psms = PsmTsvReader.ReadTsv(myFile, out warnings);  //test will fail if order of psms isn't from top to bottom row

            List<Dictionary<string, bool>> expected = new List<Dictionary<string, bool>>();
            expected.Add(new Dictionary<string, bool>() { { "b2", false }, { "b16", false }, { "y3", false }, { "y35", false } });   // no variant (7527)
            expected.Add(new Dictionary<string, bool>() { { "b2", false }, { "b3", true } });                                        // b fragments before and after variant (4221)
            expected.Add(new Dictionary<string, bool>() { { "y2", false }, { "y3", true } });                                        // y fragments before and after variant (7759)
            expected.Add(new Dictionary<string, bool>() { { "b4", true }, { "b16", true }, { "y1", false }, { "y7", false } });      // variant at 1st position (9221)
            expected.Add(new Dictionary<string, bool>() { { "b2", false }, { "b3", false }, { "y1", true }, { "y9", true } });       // variant at last position (6778)
            expected.Add(new Dictionary<string, bool>() { { "b3", false }, { "b4", true }, { "y34", false }, { "y35", true } });     // peptide comes from multiple proteins (8613)
            expected.Add(new Dictionary<string, bool>() { { "b1", true }, { "b10", true }, { "y1", true }, { "y6", true } });        // variation spans the whole peptide (8765)
            expected.Add(new Dictionary<string, bool>() { { "b2", false }, { "b5", false }, { "y1", false }, { "y30", false } });    // variant before peptide (8169)
            expected.Add(new Dictionary<string, bool>() { { "b2", false }, { "b4", false }, { "y1", false }, { "y14", false } });    // variant after peptide (6532)


            Assert.AreEqual(0, psms[0].VariantCrossingIons.Count, 
                "VariantCrossingIons should be empty for psms with no identified sequence variants");

            for (int i = 1; i < psms.Count; i++)
                foreach (MatchedFragmentIon ion in psms[i].MatchedIons)
                {
                    Assert.AreEqual(expected[i][ion.NeutralTheoreticalProduct.Annotation], psms[i].VariantCrossingIons[ion],
                       "Wrong bool value in VariantCrossingIons for ion " + ion.NeutralTheoreticalProduct.Annotation +
                       " in file " + psms.ElementAt(i).Filename + ".");
                }
                    

        }
    }
}