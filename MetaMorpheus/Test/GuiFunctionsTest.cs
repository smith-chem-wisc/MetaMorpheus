using GuiFunctions.Databases;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Net.Http;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class GuiFunctionsTest
    {

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
            if(expectedResult.Equals("\\UP000000625_reviewed.xml.gz") && isoforms) // This should only be written once, during the first test case
            {
                Console.WriteLine("Beginning Uniprot database test.");
            }
            string filename = DownloadUniProtDatabaseFunctions.GetUniprotFilename(proteomeID, reviewed, isoforms, xmlFormat, compressed);

            if (expectedResult.Equals("\\UP000000625_withUnreviewed.fasta")) // This should only be written once, during the first test case
            {
                Console.WriteLine("Finished with Uniprot database test.");
            }

            Assert.That(expectedResult, Is.EqualTo(filename));
        }

        // Occasionally the downloaded files will change, and thus the expected result will need to be updated.
        // To verify, download the database and count the number of entries.
        // The expected result will be double the number of entries, due to decoys - 9/1/23 NB
        [Test]
        [TestCase("UP000000280", true, true, true, true, "1.fasta.gz", 50)]
        [TestCase("UP000000280", true, true, true, false, "2.fasta", 50)]
        [TestCase("UP000000280", true, true, false, true, "3.fasta.gz", 50)]
        [TestCase("UP000000280", true, false, true, true, "4.xml.gz", 50)]
        [TestCase("UP000000280", false, true, true, true, "5.fasta.gz", 2)]
        [TestCase("UP000000280", true, true, false, false, "6.fasta", 50)]
        [TestCase("UP000000280", true, false, true, false, "7.xml", 50)]
        [TestCase("UP000000280", false, true, true, false, "8.fasta", 2)]
        [TestCase("UP000000280", true, false, false, false, "9.fasta", 50)]
        [TestCase("UP000000280", false, false, true, false, "10.xml", 2)]
        [TestCase("UP000000280", false, false, false, false, "11.fasta", 2)]
        public static async Task UniprotHtmlQueryTest(string proteomeID, bool reviewed, bool isoforms, bool xmlFormat, bool compressed,
           string testName, int listCount)
        {
            if (testName.Equals("1.fasta.gz")) // This should only be written once, during the first test case
            {
                Console.WriteLine("Beginning Uniprot HTML query test.");
            }

            var proteomeURL = DownloadUniProtDatabaseFunctions.GetUniProtHtmlQueryString(proteomeID, reviewed,
                isoforms, xmlFormat, compressed);

            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, $@"DatabaseTests\{testName}");

            List<Protein> reader = new List<Protein>();



            HttpClientHandler handler = new HttpClientHandler(); // without this, the download is very slow
            handler.Proxy = null;
            handler.UseProxy = false;

            var client = new HttpClient(handler); // client for using the REST Api
            var response = await client.GetAsync(proteomeURL);

            using (var file = File.Create(filePath)) // saves the file 
            {
                var content = await response.Content.ReadAsStreamAsync();
                await content.CopyToAsync(file);
            }


            if (xmlFormat)
            {
                reader = ProteinDbLoader.LoadProteinXML(proteinDbLocation: filePath, generateTargets: true, decoyType: DecoyType.Reverse,
                    allKnownModifications: null, isContaminant: false, modTypesToExclude: null, out var unknownMod);
            }
            else
            {
                reader = ProteinDbLoader.LoadProteinFasta(filePath, generateTargets: true, decoyType: DecoyType.Reverse,
                    isContaminant: false, out var unknownMod);
            }

            File.Delete(filePath);

            if(testName.Equals("11.fasta")) // only triggers for last test cases
            {
                Console.WriteLine("Finished with Uniprot HTML query test.");
            }

            Assert.That(reader.Count == listCount);
        }
    }
}