using EngineLayer;
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
        public static void TestFindVariantCrossingIons()
        {
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\VariantCrossTest.psmtsv");
            List<string> warnings = new List<string>();
            List<PsmFromTsv> psms;

            psms = PsmTsvReader.ReadTsv(myFile, out warnings);  // test will fail if the order of psms is changed to something other than top to bottom row

            // check that variant psm properties are being parsed correctly
            Assert.AreEqual("", psms[0].IdentifiedSequenceVariations);
            Assert.AreEqual("A147T", psms[1].IdentifiedSequenceVariations);

            Assert.AreEqual("", psms[0].IntersectingSequenceVariations);
            Assert.AreEqual("A147T", psms[1].IntersectingSequenceVariations);

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
                       "VariantCrossingIons should contain ion " + expectedIon + " in file " + psms[i].Filename + ".");
                foreach (string actualIon in actualIons)
                    Assert.IsTrue(expected[i].Contains(actualIon),
                        "VariantCrossingIons should not contain ion " + actualIon + " in file " + psms[i].Filename + ".");
            }
        }
    }
}