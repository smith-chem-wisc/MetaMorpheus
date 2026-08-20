using EngineLayer;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer.DatabaseLoading;
using TaskLayer;
using UsefulProteomicsDatabases;
using Transcriptomics.Digestion;

namespace Test.DatabaseTests
{
    [TestFixture]
    public class ProteinLoaderTest
    {
        [Test]
        public void ReadEmptyFasta()
        {
            new ProteinLoaderTask("").Run(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "empty.fa"));
        }

        [Test]
        public void ReadFastaWithEmptyEntry()
        {
            new ProteinLoaderTask("").Run(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "oneEmptyEntry.fa"));
        }

        [Test]
        public void TestProteinLoad()
        {
            new ProteinLoaderTask("").Run(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fasta"));
            new ProteinLoaderTask("").Run(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fa"));
            new ProteinLoaderTask("").Run(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fasta.gz"));
            new ProteinLoaderTask("").Run(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fa.gz"));
        }

        [Test]
        public void WriteTargetDecoyFasta_WhenEnabled_CreatesFastaFile()
        {
            // Arrange
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestTargetDecoyOutput");
            Directory.CreateDirectory(outputFolder);

            string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fasta");
            var dbForTask = new List<DbForTask> { new DbForTask(dbPath, false) };
            var commonParameters = new CommonParameters();

            // Act
            var loader = new DatabaseLoadingEngine(
                commonParameters,
                [],
                [],
                dbForTask,
                "TestTask",
                DecoyType.Reverse,
                true,
                null,
                TargetContaminantAmbiguity.RemoveContaminant,
                writeTargetDecoyFasta: true,
                outputFolder: outputFolder
            );
            var results = (DatabaseLoadingEngineResults)loader.Run()!;

            // Assert
            string fastaPath = Path.Combine(outputFolder, "TargetDecoy.fasta");
            Assert.That(File.Exists(fastaPath), Is.True, "TargetDecoy.fasta should be created");
            var fileContent = File.ReadAllText(fastaPath);
            Assert.That(fileContent, Is.Not.Empty, "FASTA file should not be empty");

            // Cleanup
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public void WriteTargetDecoyFasta_WhenEnabled_CreatesFastaFile_RNA()
        {
            // Arrange
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "TestTargetDecoyOutput");
            Directory.CreateDirectory(outputFolder);

            string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "20mer1.fasta");
            var dbForTask = new List<DbForTask> { new DbForTask(dbPath, false) };
            var commonParameters = new CommonParameters(digestionParams: new RnaDigestionParams());


            // Act
            var loader = new DatabaseLoadingEngine(
                commonParameters,
                [],
                [],
                dbForTask,
                "TestTask",
                DecoyType.Reverse,
                true,
                null,
                TargetContaminantAmbiguity.RemoveContaminant,
                writeTargetDecoyFasta: true,
                outputFolder: outputFolder
            );
            var results = (DatabaseLoadingEngineResults)loader.Run()!;

            // Assert
            string fastaPath = Path.Combine(outputFolder, "TargetDecoy.fasta");
            Assert.That(File.Exists(fastaPath), Is.True, "TargetDecoy.fasta should be created");
            var fileContent = File.ReadAllText(fastaPath);
            Assert.That(fileContent, Is.Not.Empty, "FASTA file should not be empty");

            // Cleanup
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public void WriteTargetDecoyFasta_WhenDisabled_DoesNotCreateFile()
        {
            // Arrange
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestNoOutput");
            Directory.CreateDirectory(outputFolder);

            string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fasta");
            var dbForTask = new List<DbForTask> { new DbForTask(dbPath, false) };

            // Act
            var loader = new DatabaseLoadingEngine(
                new CommonParameters(),
                [],
                [],
                dbForTask,
                "TestTask",
                DecoyType.Reverse,
                true,
                null,
                TargetContaminantAmbiguity.RemoveContaminant,
                writeTargetDecoyFasta: false,
                outputFolder: outputFolder
            );
            loader.Run();

            // Assert
            string fastaPath = Path.Combine(outputFolder, "TargetDecoy.fasta");
            Assert.That(File.Exists(fastaPath), Is.False, "TargetDecoy.fasta should NOT be created when disabled");

            // Cleanup
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public void WriteTargetDecoyFasta_WithNullOutputFolder_DoesNotThrow()
        {
            // Arrange
            string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fasta");
            var dbForTask = new List<DbForTask> { new DbForTask(dbPath, false) };

            // Act & Assert - Should not throw
            Assert.DoesNotThrow(() =>
            {
                var loader = new DatabaseLoadingEngine(
                    new CommonParameters(),
                    [],
                    [],
                    dbForTask,
                    "TestTask",
                    DecoyType.Reverse,
                    true,
                    null,
                    TargetContaminantAmbiguity.RemoveContaminant,
                    writeTargetDecoyFasta: true,
                    outputFolder: null
                );
                loader.Run();
            });
        }

        [Test]
        public void WriteTargetDecoyFasta_ContainsBothTargetsAndDecoys()
        {
            // Arrange
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestTargetDecoyContent");
            Directory.CreateDirectory(outputFolder);

            string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fasta");
            var dbForTask = new List<DbForTask> { new DbForTask(dbPath, false) };

            // Act
            var loader = new DatabaseLoadingEngine(
                new CommonParameters(),
                [],
                [],
                dbForTask,
                "TestTask",
                DecoyType.Reverse,
                true,
                null,
                TargetContaminantAmbiguity.RemoveContaminant,
                writeTargetDecoyFasta: true,
                outputFolder: outputFolder
            );
            var results = (DatabaseLoadingEngineResults)loader.Run()!;

            // Assert
            string fastaPath = Path.Combine(outputFolder, "TargetDecoy.fasta");
            var lines = File.ReadAllLines(fastaPath);

            // Count headers (lines starting with >)
            int headerCount = lines.Count(l => l.StartsWith(">"));
            Assert.That(headerCount, Is.EqualTo(results.BioPolymers.Count),
                "FASTA should contain all loaded biopolymers");

            // Verify at least some decoys exist (lines containing "DECOY" or similar pattern)
            bool hasDecoys = lines.Any(l => l.Contains("DECOY"));
            Assert.That(hasDecoys, Is.True, "FASTA should contain decoy sequences");

            // Cleanup
            Directory.Delete(outputFolder, true);
        }

       

        [Test]
        public static void CatchesError()
        {
            string badOutPath = @"Z:\This\Path\Does\Not\Exist";
            string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fasta");
            var dbForTask = new List<DbForTask> { new DbForTask(dbPath, false) };


            var loader = new DatabaseLoadingEngine(
                new CommonParameters(),
                [],
                [],
                dbForTask,
                "TestTask",
                DecoyType.Reverse,
                true,
                null,
                TargetContaminantAmbiguity.RemoveContaminant,
                writeTargetDecoyFasta: true,
                outputFolder: badOutPath
            );

            try
            {
                var results = (DatabaseLoadingEngineResults)loader.Run()!;
            }
            catch(System.Exception ex)
            {
                Assert.Fail("ProteinLoaderTask threw an exception: " + ex.Message);
            }
        }

        /// <summary>
        /// Writes a minimal protein database with one "modified residue" feature carrying the given description.
        /// Deliberately omits the &lt;modification&gt; blocks a MetaMorpheus-written database embeds, because those
        /// are loaded by GetPtmListFromProteinXml and would make any description resolvable.
        /// </summary>
        private static string WriteDbWithModificationDescription(string fileName, string modificationDescription)
        {
            return WriteDbWithModificationDescriptions(fileName, modificationDescription);
        }

        /// <summary>
        /// As above, but one feature per description, all annotated on the same serine, so that a single database
        /// can carry more unmatched modifications than the warning names individually.
        /// </summary>
        private static string WriteDbWithModificationDescriptions(string fileName, params string[] modificationDescriptions)
        {
            string features = string.Concat(modificationDescriptions.Select(description =>
                $"    <feature type=\"modified residue\" description=\"{description}\">\n" +
                "      <location>\n" +
                "        <position position=\"3\" />\n" +
                "      </location>\n" +
                "    </feature>\n"));

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, fileName);
            File.WriteAllText(path,
                "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n" +
                "<mzLibProteinDb>\n" +
                "  <entry>\n" +
                "  <accession>P12345</accession>\n" +
                "    <name>TEST_PROTEIN</name>\n" +
                "    <protein><recommendedName><fullName>Test protein</fullName></recommendedName></protein>\n" +
                features +
                "    <sequence length=\"11\">MPSPEPTIDEK</sequence>\n" +
                "  </entry>\n" +
                "</mzLibProteinDb>\n");
            return path;
        }

        /// <summary>
        /// A modification annotated in a database but absent from the known modifications is skipped by mzLib, which
        /// records it in its unknownModifications dictionary. That dictionary used to be discarded here, so the
        /// annotation vanished with nothing said about it (mzLib #417). It must now reach the user as a warning.
        /// </summary>
        [Test]
        public static void UnmatchedModificationInDatabaseIsWarnedAbout()
        {
            string dbPath = WriteDbWithModificationDescription("unknownModDb.xml", "Nonexistent modification on S");

            DatabaseLoadingEngine.LoadBioPolymers("TestTask", new List<DbForTask> { new DbForTask(dbPath, false) },
                true, DecoyType.None, new List<string> { "UniProt" }, new CommonParameters(), out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(errors.Any(p => p.Contains("Nonexistent modification on S")),
                    $"expected the unmatched modification to be named; got: {string.Join(" | ", errors)}");
                Assert.That(errors.Any(p => p.Contains("unknownModDb.xml")),
                    "expected the warning to name the database it came from");
            });
        }

        /// <summary>
        /// The counterpart: a modification that does resolve must not be reported. Without this, the warning above
        /// could be satisfied by warning unconditionally, which would be worse than the silence it replaces.
        /// </summary>
        [Test]
        public static void ResolvableModificationProducesNoUnmatchedModificationWarning()
        {
            string dbPath = WriteDbWithModificationDescription("knownModDb.xml", "Phosphoserine on S");

            DatabaseLoadingEngine.LoadBioPolymers("TestTask", new List<DbForTask> { new DbForTask(dbPath, false) },
                true, DecoyType.None, new List<string> { "UniProt" }, new CommonParameters(), out var errors);

            Assert.That(errors.Any(p => p.Contains("could not be matched")), Is.False,
                $"no unmatched-modification warning expected; got: {string.Join(" | ", errors)}");
        }

        /// <summary>
        /// A database annotating more unmatched modifications than the warning names individually: the first five
        /// in ordinal order are named and the remainder is reported as a count, so a database using many unknown
        /// ids produces one bounded message rather than an unbounded list.
        /// </summary>
        [Test]
        public static void UnmatchedModificationWarningNamesFiveThenCountsTheRest()
        {
            string[] descriptions = Enumerable.Range(1, 7)
                .Select(i => $"Nonexistent modification {i} on S")
                .ToArray();
            string dbPath = WriteDbWithModificationDescriptions("manyUnknownModsDb.xml", descriptions);

            DatabaseLoadingEngine.LoadBioPolymers("TestTask", new List<DbForTask> { new DbForTask(dbPath, false) },
                true, DecoyType.None, new List<string> { "UniProt" }, new CommonParameters(), out var errors);

            string warning = errors.SingleOrDefault(p => p.Contains("could not be matched"));
            Assert.That(warning, Is.Not.Null,
                $"expected exactly one unmatched-modification warning; got: {string.Join(" | ", errors)}");
            Assert.Multiple(() =>
            {
                // the count covers every unmatched id, not just the named ones
                Assert.That(warning, Does.Contain("7 annotated modification(s)"));
                Assert.That(warning, Does.Contain("and 2 more"));
                Assert.That(warning, Does.Contain("'Nonexistent modification 1 on S'"));
                Assert.That(warning, Does.Contain("'Nonexistent modification 5 on S'"));
                Assert.That(warning, Does.Not.Contain("'Nonexistent modification 6 on S'"));
            });
        }

        public class ProteinLoaderTask : MetaMorpheusTask
        {
            public ProteinLoaderTask(string x)
                : this()
            { }

            protected ProteinLoaderTask()
                : base(MyTask.Search)
            { }

            public void Run(string dbPath)
            {
                RunSpecific("", new List<DbForTask> { new DbForTask(dbPath, false) }, null, "", null);
            }

            protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
            {
                var dbLoader = new DatabaseLoadingEngine(new(), [], [], dbFilenameList, taskId, DecoyType.None);
                var results = dbLoader.Run();
                return null;
            }
        }
    }
}