using GuiFunctions.Databases;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Text;

namespace Test
{
    [TestFixture]
    public static class GuiFunctionsTest
    {
        [Test]
        [TestCase("UP000000625", true, false, false, true, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625 reviewed:yes&compress=yes&format=fasta&include:no")]
        [TestCase("UP000000625", false, true, false, true, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625&compress=yes&format=fasta&include:yes")]
        [TestCase("UP000000625", false, true, false, false, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625&compress=no&format=fasta&include:yes")]
        [TestCase("UP000000625", false, false, false, true, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625&compress=yes&format=fasta&include:no")]
        [TestCase("UP000000625", false, false, false, false, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625&compress=no&format=fasta&include:no")]
        [TestCase("UP000000625", false, false, true, true, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625&compress=yes&format=xml&include:no")]
        [TestCase("UP000000625", false, false, true, false, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625&compress=no&format=xml&include:no")]
        [TestCase("UP000000625", true, true, true, true, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625 reviewed:yes&compress=yes&format=xml&include:yes")]
        [TestCase("UP000000625", true, true, true, false, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625 reviewed:yes&compress=no&format=xml&include:yes")]
        [TestCase("UP000000625", true, true, false, true, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625 reviewed:yes&compress=yes&format=fasta&include:yes")]
        [TestCase("UP000000625", true, false, true, true, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625 reviewed:yes&compress=yes&format=xml&include:no")]
        [TestCase("UP000000625", false, true, true, true, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625&compress=yes&format=xml&include:yes")]
        [TestCase("UP000000625", true, true, false, false, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625 reviewed:yes&compress=no&format=fasta&include:yes")]
        [TestCase("UP000000625", true, false, true, false, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625 reviewed:yes&compress=no&format=xml&include:no")]
        [TestCase("UP000000625", false, true, true, false, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625&compress=no&format=xml&include:yes")]
        [TestCase("UP000000625", true, false, false, false, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625 reviewed:yes&compress=no&format=fasta&include:no")]
        [TestCase("UP000000625", false, false, true, false, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625&compress=no&format=xml&include:no")]
        [TestCase("UP000000625", false, false, false, false, "https://www.uniprot.org/uniprot/?query=proteome:UP000000625&compress=no&format=fasta&include:no")]
        public static void TestGetUniProtHtmlQueryString(string proteomeID, bool reviewed, bool isoforms, bool xmlFormat, bool compressed, string expectedResult)
        {
            string queryString = DownloadUniProtDatabaseFunctions.GetUniProtHtmlQueryString(proteomeID, reviewed, isoforms, xmlFormat, compressed);
            Assert.AreEqual(expectedResult, queryString);
        }

        [Test]
        [TestCase("UP000000625", true, true, true, true, "\\UP000000625_reviewed_isoform.xml.gz")]
        [TestCase("UP000000625", true, true, true, false, "\\UP000000625_reviewed_isoform.xml")]
        [TestCase("UP000000625", true, true, false, true, "\\UP000000625_reviewed_isoform.fasta.gz")]
        [TestCase("UP000000625", true, false, true, true, "\\UP000000625_reviewed.xml.gz")]
        [TestCase("UP000000625", false, true, true, true, "\\UP000000625_withUnreviewed_isoform.xml.gz")]
        [TestCase("UP000000625", true, true, false, false, "\\UP000000625_reviewed_isoform.fasta")]
        [TestCase("UP000000625", true, false, true, false, "\\UP000000625_reviewed.xml")]
        [TestCase("UP000000625", false, true, true, false, "\\UP000000625_withUnreviewed_isoform.xml")]
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
