using GuiFunctions.Databases;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    public static class GuiFunctionsTest
    {
        [Test]
        [TestCase("UP000000625", true, false, false, true, "https://rest.uniprot.org/uniprotkb/stream?format=fasta&reviewed=true&compressed=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", false, true, false, true, "https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true&compressed=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", false, true, false, false, "https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", false, false, false, true, "https://rest.uniprot.org/uniprotkb/stream?format=fasta&compressed=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", false, false, false, false, "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", false, false, true, true, "https://rest.uniprot.org/uniprotkb/stream?format=xml&compressed=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", false, false, true, false, "https://rest.uniprot.org/uniprotkb/stream?format=xml&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", true, true, true, true, "https://rest.uniprot.org/uniprotkb/stream?format=xml&reviewed=true&compressed=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", true, true, true, false, "https://rest.uniprot.org/uniprotkb/stream?format=xml&reviewed=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", true, true, false, true, "https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true&reviewed=true&compressed=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", true, false, true, true, "https://rest.uniprot.org/uniprotkb/stream?format=xml&reviewed=true&compressed=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", false, true, true, true, "https://rest.uniprot.org/uniprotkb/stream?format=xml&compressed=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", true, true, false, false, "https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true&reviewed=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", true, false, true, false, "https://rest.uniprot.org/uniprotkb/stream?format=xml&reviewed=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", false, true, true, false, "https://rest.uniprot.org/uniprotkb/stream?format=xml&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", true, false, false, false, "https://rest.uniprot.org/uniprotkb/stream?format=fasta&reviewed=true&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", false, false, true, false, "https://rest.uniprot.org/uniprotkb/stream?format=xml&query=%28%28proteome%3AUP000000625%29%29")]
        [TestCase("UP000000625", false, false, false, false, "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28proteome%3AUP000000625%29%29")]
        public static void TestGetUniProtHtmlQueryString(string proteomeID, bool reviewed, bool isoforms, bool xmlFormat, bool compressed, string expectedResult)
        {
            string queryString = DownloadUniProtDatabaseFunctions.GetUniProtHtmlQueryString(proteomeID, reviewed, isoforms, xmlFormat, compressed);
            Assert.AreEqual(expectedResult, queryString);
        }

        [Test]
        [TestCase("UP000000625", true, true, true, true, "\\UP000000625_reviewed.xml.gz")]
        [TestCase("UP000000625", true, true, true, false, "\\UP000000625_reviewed.xml")]
        [TestCase("UP000000625", true, true, false, true, "\\UP000000625_reviewed_isoform.fasta.gz")]
        [TestCase("UP000000625", true, false, true, true, "\\UP000000625_reviewed.xml.gz")]
        [TestCase("UP000000625", false, true, true, true, "\\UP000000625_withUnreviewed.xml.gz")]
        [TestCase("UP000000625", true, true, false, false, "\\UP000000625_reviewed_isoform.fasta")]
        [TestCase("UP000000625", true, false, true, false, "\\UP000000625_reviewed.xml")]
        [TestCase("UP000000625", false, true, true, false, "\\UP000000625_withUnreviewed.xml")]
        [TestCase("UP000000625", true, false, false, false, "\\UP000000625_reviewed.fasta")]
        [TestCase("UP000000625", false, false, true, false, "\\UP000000625_withUnreviewed.xml")]
        [TestCase("UP000000625", false, false, false, false, "\\UP000000625_withUnreviewed.fasta")]
        public static void TestGetUniprotFilename(string proteomeID, bool reviewed, bool isoforms, bool xmlFormat, bool compressed, string expectedResult)
        {
            string filename = DownloadUniProtDatabaseFunctions.GetUniprotFilename(proteomeID, reviewed, isoforms, xmlFormat, compressed);
            Assert.AreEqual(expectedResult, filename);
        }
    }
}