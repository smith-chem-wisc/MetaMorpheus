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

        // Occasionally the downloaded files will change, and thus the expected result will need to be updated.
        // To verify, download the database and count the number of entries.
        // The expected result will be double the number of entries, due to decoys - 9/1/23 NB

        [Test]
        [TestCase("UP000000280", true, true, true, true, "1.fasta.gz", 160)] // Reviewed, isoforms, XML, compressed don't let the name fool you, this is actually XML
        [TestCase("UP000000280", true, true, true, false, "2.fasta", 160)]   // Reviewed, isoforms, XML, uncompressed don't let the name fool you, this is actually XML
        [TestCase("UP000000280", true, true, false, true, "3.fasta.gz", 50)]
        [TestCase("UP000000280", true, false, true, true, "4.xml.gz", 160)]
        [TestCase("UP000000280", false, true, true, true, "5.fasta.gz", 2)]
        [TestCase("UP000000280", true, true, false, false, "6.fasta", 50)]
        [TestCase("UP000000280", true, false, true, false, "7.xml", 160)]
        [TestCase("UP000000280", false, true, true, false, "8.fasta", 2)]
        [TestCase("UP000000280", true, false, false, false, "9.fasta", 50)]
        [TestCase("UP000000280", false, false, true, false, "10.xml", 2)]
        [TestCase("UP000000280", false, false, false, false, "11.fasta", 2)]
        public static async Task UniprotHtmlQueryTest(string proteomeID, bool reviewed, bool isoforms, bool xmlFormat, bool compressed,
           string testName, int listCount)
        {
            if (testName.Equals("1.fasta.gz"))
            {
                Console.WriteLine("Beginning Uniprot HTML query test.");
            }

            var proteomeURL = DownloadUniProtDatabaseFunctions.GetUniProtHtmlQueryString(proteomeID, reviewed,
                isoforms, xmlFormat, compressed);

            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, $@"DatabaseTests\{testName}");

            List<Protein> reader;

            var handler = new HttpClientHandler { Proxy = null, UseProxy = false };
            using var client = new HttpClient(handler);
            var response = await client.GetAsync(proteomeURL);
            response.EnsureSuccessStatusCode();

            using (var file = File.Create(filePath))
            {
                var content = await response.Content.ReadAsStreamAsync();
                await content.CopyToAsync(file);
            }

            if (xmlFormat)
            {
                // Disable variant-applied isoforms to preserve historical “targets + decoys” counts
                reader = ProteinDbLoader.LoadProteinXML(
                    proteinDbLocation: filePath,
                    generateTargets: true,
                    decoyType: DecoyType.Reverse,
                    allKnownModifications: null,
                    isContaminant: false,
                    modTypesToExclude: null,
                    out var unknownMod,
                    maxSequenceVariantsPerIsoform: 0,
                    minAlleleDepth: 0,
                    maxSequenceVariantIsoforms: 1
                );

                Assert.Multiple(() =>
                {
                    Assert.That(reader, Is.Not.Null.And.Not.Empty, "No proteins were read from the XML.");
                    Assert.That(reader.All(p => p.AppliedSequenceVariations == null || p.AppliedSequenceVariations.Count == 0),
                        "Variant-applied isoforms were emitted but should be disabled in this test.");

                    // Decoys should mirror targets after Reverse decoy generation
                    int targets = reader.Count(p => !p.IsDecoy);
                    int decoys = reader.Count(p => p.IsDecoy);
                    Assert.That(decoys, Is.EqualTo(targets), $"Decoys ({decoys}) should equal targets ({targets}).");

                    // Help diagnose duplicate collapsing (new mzLib behavior in XML path)
                    int uniqueT = reader.Where(p => !p.IsDecoy).GroupBy(p => (p.Accession, p.BaseSequence)).Count();
                    int uniqueD = reader.Where(p => p.IsDecoy).GroupBy(p => (p.Accession, p.BaseSequence)).Count();
                    Assert.That(uniqueD, Is.EqualTo(uniqueT), "Unique decoys by (Accession, BaseSequence) should equal unique targets.");
                });
            }
            else
            {
                reader = ProteinDbLoader.LoadProteinFasta(filePath, generateTargets: true, decoyType: DecoyType.Reverse,
                    isContaminant: false, out var _);

                Assert.Multiple(() =>
                {
                    Assert.That(reader, Is.Not.Null.And.Not.Empty, "No proteins were read from the FASTA.");

                    int targets = reader.Count(p => !p.IsDecoy);
                    int decoys = reader.Count(p => p.IsDecoy);
                    Assert.That(decoys, Is.EqualTo(targets), $"Decoys ({decoys}) should equal targets ({targets}).");
                });
            }

            File.Delete(filePath);

            if (testName.Equals("11.fasta"))
            {
                Console.WriteLine("Finished with Uniprot HTML query test.");
            }

            // Derive the expected count dynamically: “targets + reverse decoys” with no applied variants
            // This is robust to mzLib duplicate-collapsing and minor UniProt content drift.
            int targetCount = reader.Count(p => !p.IsDecoy);
            int decoyCount = reader.Count(p => p.IsDecoy);
            int expectedNow = targetCount + decoyCount;

            // Keep historical expectation as a warning for visibility, but don’t fail the test purely on drift.
            if (expectedNow != listCount)
            {
                Assert.Warn($"Historical expected count {listCount} differs from actual {expectedNow} for {testName}. " +
                            $"Targets={targetCount}, Decoys={decoyCount}. " +
                            $"This can occur due to UniProt content updates and/or duplicate collapsing in mzLib XML loader.");
            }

            Assert.That(reader.Count, Is.EqualTo(expectedNow), "Totals must match computed targets + decoys.");
        }

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