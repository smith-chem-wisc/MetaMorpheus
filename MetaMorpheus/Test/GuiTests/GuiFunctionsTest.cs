using GuiFunctions.Databases;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Net.Http;
using System.Threading.Tasks;
using GuiFunctions;
using UsefulProteomicsDatabases;
using System.Linq;

namespace Test.GuiTests
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
            if (expectedResult.Equals("\\UP000000625_reviewed.xml.gz") && isoforms) // This should only be written once, during the first test case
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

        // Verifies GetUniProtHtmlQueryString builds the exact REST URL the GUI submits to
        // rest.uniprot.org for every reviewed/isoform/format/compressed combination.
        // This is the offline replacement for the assertions that used to depend on a live
        // UniProt download; the live download itself is exercised by UniProtDownloadCanary below.
        [Test]
        [TestCase("UP000001207", true, true, true, true, "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=xml&query=reviewed:true AND proteome:UP000001207")]
        [TestCase("UP000001207", true, true, true, false, "https://rest.uniprot.org/uniprotkb/stream?format=xml&query=reviewed:true AND proteome:UP000001207")]
        [TestCase("UP000001207", true, true, false, true, "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&includeIsoform=true&query=reviewed:true AND proteome:UP000001207")]
        [TestCase("UP000001207", true, false, true, true, "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=xml&query=reviewed:true AND proteome:UP000001207")]
        [TestCase("UP000001207", false, true, true, true, "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=xml&query=proteome:UP000001207")]
        [TestCase("UP000001207", true, true, false, false, "https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true&query=reviewed:true AND proteome:UP000001207")]
        [TestCase("UP000001207", true, false, true, false, "https://rest.uniprot.org/uniprotkb/stream?format=xml&query=reviewed:true AND proteome:UP000001207")]
        [TestCase("UP000001207", false, true, true, false, "https://rest.uniprot.org/uniprotkb/stream?format=xml&query=proteome:UP000001207")]
        [TestCase("UP000001207", true, false, false, false, "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=reviewed:true AND proteome:UP000001207")]
        [TestCase("UP000001207", false, false, true, false, "https://rest.uniprot.org/uniprotkb/stream?format=xml&query=proteome:UP000001207")]
        [TestCase("UP000001207", false, false, false, false, "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=proteome:UP000001207")]
        public static void TestGetUniProtHtmlQueryString(string proteomeID, bool reviewed, bool isoforms, bool xmlFormat, bool compressed, string expectedUrl)
        {
            string url = DownloadUniProtDatabaseFunctions.GetUniProtHtmlQueryString(proteomeID, reviewed, isoforms, xmlFormat, compressed);
            Assert.That(url, Is.EqualTo(expectedUrl));
        }

        // Offline "local mock" of a downloaded UniProt proteome. Rather than streaming from
        // rest.uniprot.org (which is unreliable in CI), we load committed UniProtKB-format
        // fixtures through the same ProteinDbLoader path the download uses and confirm that
        // targets + reverse decoys are generated (count == 2 * entries). The fixtures live in
        // DatabaseTests\UniProtMock and contain 3 entries each.
        [Test]
        [TestCase(@"DatabaseTests\UniProtMock\phi29_mock_reviewed.xml", true, 6)]   // 3 entries * 2 (target + decoy)
        [TestCase(@"DatabaseTests\UniProtMock\phi29_mock_reviewed.fasta", false, 6)] // 3 entries * 2 (target + decoy)
        public static void LoadUniProtProteomeFromLocalFile_GeneratesTargetsAndDecoys(string relativePath, bool xmlFormat, int expectedCount)
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, relativePath);

            List<Protein> reader;
            if (xmlFormat)
            {
                reader = ProteinDbLoader.LoadProteinXML(proteinDbLocation: filePath, generateTargets: true, decoyType: DecoyType.Reverse,
                    allKnownModifications: null, isContaminant: false, modTypesToExclude: null, out _, maxHeterozygousVariants: 0);
            }
            else
            {
                reader = ProteinDbLoader.LoadProteinFasta(filePath, generateTargets: true, decoyType: DecoyType.Reverse,
                    isContaminant: false, out _);
            }

            Assert.That(reader.Count, Is.EqualTo(expectedCount));
        }

        // Live canary for the UniProt REST download used by the "Download UniProt Database" GUI.
        // Tagged [Category("ExternalService")] (+ "UniProt") so it runs in the dedicated, non-blocking
        // external-service CI job rather than the required unit-test run. It confirms rest.uniprot.org
        // is reachable, that GetUniProtHtmlQueryString still yields a working URL, and that
        // ProteinDbLoader can parse the response. Via ExternalServiceTestHelper.RunAsync, a UniProt
        // outage (transport error, 5xx, or the HTTP-200 "Error encountered when streaming data" body)
        // is reported as Skipped rather than Failed - "we tried, UniProt was down, don't worry" - while
        // a genuine contract break (nothing parses, URL rejected) still fails. The assertion is
        // intentionally loose (proteins parsed > 0) so a UniProt content revision is not a red herring.
        // UP000001207 = Bacillus phage phi29 (small reference proteome).
        [Test]
        [Category("ExternalService")]
        [Category("UniProt")]
        [TestCase("UP000001207", true, false, true, false)]   // reviewed, XML
        [TestCase("UP000001207", false, false, false, false)] // all, FASTA
        public static Task UniProtDownloadCanary(string proteomeID, bool reviewed, bool isoforms, bool xmlFormat, bool compressed) =>
            ExternalServiceTestHelper.RunAsync("UniProt", async () =>
            {
                var proteomeURL = DownloadUniProtDatabaseFunctions.GetUniProtHtmlQueryString(proteomeID, reviewed,
                    isoforms, xmlFormat, compressed);

                var extension = (xmlFormat ? ".xml" : ".fasta") + (compressed ? ".gz" : "");
                var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, $@"DatabaseTests\uniprot_canary{extension}");

                HttpClientHandler handler = new HttpClientHandler(); // without this, the download is very slow
                handler.Proxy = null;
                handler.UseProxy = false;

                using var client = new HttpClient(handler); // client for using the REST Api
                var response = await client.GetAsync(proteomeURL);

                var bytes = await response.Content.ReadAsByteArrayAsync();
                // Sniff a text prefix (harmless on gzipped payloads) so a service outage is skipped, not failed.
                var textPreview = System.Text.Encoding.UTF8.GetString(bytes, 0, Math.Min(bytes.Length, 512));
                ExternalServiceTestHelper.ThrowIfUnavailable(response, textPreview);

                File.WriteAllBytes(filePath, bytes); // saves the file

                List<Protein> reader;
                if (xmlFormat)
                {
                    reader = ProteinDbLoader.LoadProteinXML(proteinDbLocation: filePath, generateTargets: true, decoyType: DecoyType.Reverse,
                        allKnownModifications: null, isContaminant: false, modTypesToExclude: null, out _, maxHeterozygousVariants: 0);
                }
                else
                {
                    reader = ProteinDbLoader.LoadProteinFasta(filePath, generateTargets: true, decoyType: DecoyType.Reverse,
                        isContaminant: false, out _);
                }

                File.Delete(filePath);

                Assert.That(reader.Count, Is.GreaterThan(0),
                    "UniProt returned no parseable proteins - the REST endpoint may be down or the query URL/format has changed.");
            });

        [Test]
        public static void TestFileLoadingWithDuplicateFiles()
        {
            var metaDrawLogic = new MetaDrawLogic();
            metaDrawLogic.SpectraFilePaths.Add(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\mouseOne.mzML"));
            metaDrawLogic.SpectraFilePaths.Add(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\mouseOne.mzML"));

            var warnings = metaDrawLogic.LoadFiles(true, false);
            Assert.That(warnings.Count, Is.EqualTo(0));
            Assert.That(metaDrawLogic.MsDataFiles.Count, Is.EqualTo(1));
        }

        [Test]
        public static void TestFileLoadingTimsTof()
        {
            var metaDrawLogic = new MetaDrawLogic();
            metaDrawLogic.SpectraFilePaths.Add(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\snippet.d"));

            var warnings = metaDrawLogic.LoadFiles(true, false);
            Assert.That(warnings.Count, Is.EqualTo(0));
            Assert.That(metaDrawLogic.MsDataFiles.Count, Is.EqualTo(1));
            Assert.That(metaDrawLogic.MsDataFiles.First().Value.Scans.Length, Is.GreaterThan(0)); // Check that scans have already been read in
        }
    }
}
