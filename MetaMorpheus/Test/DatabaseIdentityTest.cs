using EngineLayer.DatabaseLoading;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Security.Cryptography;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Every task's outputs say exactly which database files it searched: their content (SHA-256 and size), not
    /// only their names, and for a GPTMD-written database, the databases it was written from.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class DatabaseIdentityTest
    {
        private static string Sha256Of(string path) =>
            Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path))).ToLowerInvariant();

        [Test]
        public static void RecordIdentity_IsTheSha256AndSizeOfTheFileAsGiven()
        {
            string database = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mouseOne.xml");
            var db = new DbForTask(database, false);
            Assert.That(db.Sha256, Is.Null, "nothing is known before the database is loaded");

            db.RecordIdentity();

            Assert.That(db.Sha256, Is.EqualTo(Sha256Of(database)));
            Assert.That(db.SizeBytes, Is.EqualTo(new FileInfo(database).Length));
            Assert.That(db.IdentityError, Is.Null);
            Assert.That(db.IdentityText(), Is.EqualTo($"SHA-256 {Sha256Of(database)}, {new FileInfo(database).Length} bytes"));
            Assert.That(db.DerivedFromText(), Is.Null, "a database the user supplied was not derived from another");
        }

        /// <summary>The same path holding new content is a different database, and the identity follows the content.</summary>
        [Test]
        public static void RecordIdentity_FollowsTheContentNotThePath()
        {
            string dir = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseIdentity_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(dir);
            try
            {
                string path = Path.Combine(dir, "proteins.fasta");
                File.WriteAllText(path, ">sp|P1|ONE\nPEPTIDEK\n");
                var db = new DbForTask(path, false);
                db.RecordIdentity();
                string first = db.Sha256;

                File.WriteAllText(path, ">sp|P2|TWO\nKEDITPEP\n");
                db.RecordIdentity();

                Assert.That(db.Sha256, Is.Not.EqualTo(first));
                Assert.That(db.Sha256, Is.EqualTo(Sha256Of(path)));
            }
            finally
            {
                Directory.Delete(dir, true);
            }
        }

        /// <summary>A file that cannot be read is recorded as unidentified, with the reason, instead of throwing.</summary>
        [Test]
        public static void RecordIdentity_MissingFile_RecordsWhyAndDoesNotThrow()
        {
            var db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, "no_such_database.xml"), false);

            Assert.DoesNotThrow(db.RecordIdentity);

            Assert.That(db.Sha256, Is.Null);
            Assert.That(db.SizeBytes, Is.Null);
            Assert.That(db.IdentityError, Is.Not.Null.And.Not.Empty);
            Assert.That(db.IdentityText(), Does.StartWith("SHA-256 unavailable ("));
        }

        /// <summary>
        /// GPTMD, then a search of the database GPTMD wrote. The search's results.txt and prose must name that
        /// written file's content, and say which database it was written from and that one's content.
        /// </summary>
        [Test]
        public static void GptmdThenSearch_OutputsNameEachDatabaseByContentAndParent()
        {
            string spectra = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mouseOne.mzML");
            string database = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mouseOne.xml");
            string output = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseIdentity_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(output);
            try
            {
                var tasks = new List<(string, MetaMorpheusTask)>
                {
                    ("gptmd", new GptmdTask()),
                    ("search", new SearchTask { SearchParameters = new SearchParameters { DoLabelFreeQuantification = false } }),
                };
                new EverythingRunnerEngine(tasks, new List<string> { spectra }, new List<DbForTask> { new DbForTask(database, false) }, output).Run();

                string parentHash = Sha256Of(database);
                string gptmdDb = Directory.GetFiles(Path.Combine(output, "gptmd"), "*GPTMD.xml").Single();
                string gptmdHash = Sha256Of(gptmdDb);

                // GPTMD's own results name the database it wrote and where it came from.
                string gptmdResults = File.ReadAllText(Path.Combine(output, "gptmd", "results.txt"));
                Assert.That(gptmdResults, Does.Contain($"\t{Path.GetFileName(database)}: "));
                Assert.That(gptmdResults, Does.Contain($"\t\tPath: {database}"));
                Assert.That(gptmdResults, Does.Contain($"\t\tSHA-256 {parentHash}, "));
                Assert.That(gptmdResults, Does.Contain($"\t{Path.GetFileName(gptmdDb)}: SHA-256 {gptmdHash}, "));

                // The search loaded the written file, and says so by content and by parent.
                string[] searchResults = File.ReadAllLines(Path.Combine(output, "search", "results.txt"));
                int at = Array.FindIndex(searchResults, l => l.StartsWith($"\t{Path.GetFileName(gptmdDb)}: ", StringComparison.Ordinal));
                Assert.That(at, Is.GreaterThanOrEqualTo(0), "the search lists the GPTMD database it loaded");
                Assert.That(searchResults[at + 1], Does.Contain("Target"), "the existing counts line is unchanged and still follows the name");
                Assert.That(searchResults[at + 2], Is.EqualTo($"\t\tPath: {gptmdDb}"));
                Assert.That(searchResults[at + 3], Does.StartWith($"\t\tSHA-256 {gptmdHash}, "));
                Assert.That(searchResults[at + 4], Does.StartWith($"\t\tDerived from: {Path.GetFileName(database)} (SHA-256 {parentHash}, "));

                // The prose gives one line per database, with the same identity.
                string[] prose = File.ReadAllLines(Path.Combine(output, "search", "AutoGeneratedManuscriptProse.txt"));
                string databaseLine = prose.SkipWhile(l => l != "Databases:").Skip(1).First();
                Assert.That(databaseLine, Does.StartWith("\t" + gptmdDb + " Downloaded on: "));
                Assert.That(databaseLine, Does.Contain($"(SHA-256 {gptmdHash}, "));
                Assert.That(databaseLine, Does.Contain($"Derived from: {Path.GetFileName(database)} (SHA-256 {parentHash}, "));
            }
            finally
            {
                Directory.Delete(output, true);
            }
        }
    }
}
