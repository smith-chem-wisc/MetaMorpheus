using EngineLayer;
using GuiFunctions;
using NUnit.Framework;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class PsvTsvTest
    {
        [Test]
        public static void ReadOGlycoSinglePsms()
        {
            string psmFile = @"TestData\oglycoSinglePsms.psmtsv";
            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            Assert.AreEqual(2, parsedPsms.Count);
        }

        [Test]
        public static void ReadOGlycoPsms()
        {
            string psmFile = @"TestData\oglyco.psmtsv";
            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            Assert.AreEqual(9, parsedPsms.Count);
        }

        [Test]
        public static void ReadExcelEditedPsms()
        {
            string psmFile = @"TestData\ExcelEditedPeptide.psmtsv";
            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            Assert.AreEqual(1, parsedPsms.Count);
            IEnumerable<string> expectedIons = new string[] { "y3+1", "y4+1", "b4+1", "b5+1", "b6+1", "b8+1" };
            Assert.That(6 == parsedPsms[0].MatchedIons.Select(p => p.Annotation).Intersect(expectedIons).Count());
        }

        [Test]
        public static void MetaDrawLogicTestOglyco()
        {
            var metadrawLogic = new MetaDrawLogic();
            string psmFile = @"TestData\oglyco.psmtsv";
            metadrawLogic.PsmResultFilePaths.Add(psmFile);
            var errors = metadrawLogic.LoadFiles(false, true);

            Assert.That(!errors.Any());
        }
    }
}