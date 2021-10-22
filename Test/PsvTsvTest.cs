using EngineLayer;
using NUnit.Framework;
using System.Collections.Generic;

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
    }
}