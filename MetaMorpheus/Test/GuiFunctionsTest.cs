using GuiFunctions.Databases;
using NUnit.Framework;
using MetaMorpheusGUI.Properties;
using System.Net.Http;
using System;
using System.Drawing.Printing;
using System.Threading.Tasks;
using Microsoft.ML.Data;
using System.IO;
using System.Net;
using System.Security.Cryptography;
using System.Windows.Documents;
using System.Windows.Media.Animation;
using pepXML.Generated;
using UsefulProteomicsDatabases;

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

        [Test]
        //[TestCase("UP000325672", true, true, false, true, "1")]
        //[TestCase("UP000325672", true, true, false, false, "2")]
        //[TestCase("UP000325672", true, false, true, true, "3")]
        //[TestCase("UP000325672", true, false, true, false, "4")]
        [Parallelizable(ParallelScope.All)]
        [TestCase("UP000000280", true, true, true, true, "1.fasta.gz")]
        [TestCase("UP000000280", true, true, true, false, "2.fasta")]
        [TestCase("UP000000280", true, true, false, true, "3.fasta.gz")]
        [TestCase("UP000000280", true, false, true, true, "4.xml.gz")]
        [TestCase("UP000000280", false, true, true, true, "5.fasta.gz")]
        [TestCase("UP000000280", true, true, false, false, "6.fasta")]
        [TestCase("UP000000280", true, false, true, false, "7.xml")]
        [TestCase("UP000000280", false, true, true, false, "8.fasta")]
        [TestCase("UP000000280", true, false, false, false, "9.fasta")]
        [TestCase("UP000000280", false, false, true, false, "10.xml")]
        [TestCase("UP000000280", false, false, false, false, "11.fasta")]
        public static void UniprotHtmlQueryTest(string proteomeID, bool reviewed, bool isoforms, bool xmlFormat, bool compressed, 
            string testName)
        {
            var proteomeURL = DownloadUniProtDatabaseFunctions.GetUniProtHtmlQueryString(proteomeID, reviewed,
                isoforms, xmlFormat, compressed);

            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, $@"DatabaseTests\{testName}");

            try
            {
                using (var client = new WebClient())
                {
                    ServicePointManager.DefaultConnectionLimit = 15;
                    client.Proxy = null;
                    client.DownloadFile(proteomeURL, filePath);
                }

                if (xmlFormat)
                {
                    ProteinDbLoader.LoadProteinXML(proteinDbLocation: filePath, generateTargets: true, decoyType: DecoyType.Reverse,
                        allKnownModifications: null, isContaminant: false, modTypesToExclude: null, out var unknownMod);
                }
                else 
                {
                    ProteinDbLoader.LoadProteinFasta(filePath, generateTargets: true, decoyType: DecoyType.Reverse,
                        isContaminant: false, out var unknownMod);
                }

                File.Delete(filePath);

            }
            catch (Exception ex)
            {
                Console.WriteLine(ex);
                Console.WriteLine(proteomeURL);
                Assert.Fail();
            }
            Assert.Pass();
            
        }

        [Test]
        //[TestCase("UP000325672", true, true, false, true, "1")]
        //[TestCase("UP000325672", true, true, false, false, "2")]
        //[TestCase("UP000325672", true, false, true, true, "3")]
        //[TestCase("UP000325672", true, false, true, false, "4")]
        [Parallelizable(ParallelScope.All)]
        [TestCase("UP000000280", true, true, true, true, "1.fasta.gz")]
        [TestCase("UP000000280", true, true, true, false, "2.fasta")]
        [TestCase("UP000000280", true, true, false, true, "3.fasta.gz")]
        [TestCase("UP000000280", true, false, true, true, "4.xml.gz")]
        [TestCase("UP000000280", false, true, true, true, "5.fasta.gz")]
        [TestCase("UP000000280", true, true, false, false, "6.fasta")]
        [TestCase("UP000000280", true, false, true, false, "7.xml")]
        [TestCase("UP000000280", false, true, true, false, "8.fasta")]
        [TestCase("UP000000280", true, false, false, false, "9.fasta")]
        [TestCase("UP000000280", false, false, true, false, "10.xml")]
        [TestCase("UP000000280", false, false, false, false, "11.fasta")]
        public static async Task  UniprotHtmlQueryTest2(string proteomeID, bool reviewed, bool isoforms, bool xmlFormat, bool compressed,
           string testName)
        {
            var proteomeURL = DownloadUniProtDatabaseFunctions.GetUniProtHtmlQueryString(proteomeID, reviewed,
                isoforms, xmlFormat, compressed);

            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, $@"DatabaseTests\{testName}");

            try
            {
                HttpClientHandler handler = new HttpClientHandler();
                handler.Proxy = null;
                handler.UseProxy = false;

                var client = new HttpClient(handler);
                var response = await client.GetAsync(proteomeURL);
                using (var file = File.Create(filePath))
                {
                    var content = await response.Content.ReadAsStreamAsync();
                    await content.CopyToAsync(file);
                }

                if (xmlFormat)
                {
                    ProteinDbLoader.LoadProteinXML(proteinDbLocation: filePath, generateTargets: true, decoyType: DecoyType.Reverse,
                        allKnownModifications: null, isContaminant: false, modTypesToExclude: null, out var unknownMod);
                }
                else
                {
                    ProteinDbLoader.LoadProteinFasta(filePath, generateTargets: true, decoyType: DecoyType.Reverse,
                        isContaminant: false, out var unknownMod);
                }

                //File.Delete(filePath);

            }
            catch (Exception ex)
            {
                Console.WriteLine(ex);
                Console.WriteLine(proteomeURL);
                Assert.Fail();
            }
            Assert.Pass();

        }
    }
}