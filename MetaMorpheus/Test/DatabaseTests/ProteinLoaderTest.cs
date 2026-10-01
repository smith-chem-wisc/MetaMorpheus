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

        private AnalyteType _analyteTypeBeforeTest;

        private static string ScratchDirectory => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "ProteinLoaderTestScratch", TestContext.CurrentContext.Test.MethodName);

        [SetUp]
        public void RecordProcessWideState() => _analyteTypeBeforeTest = GlobalVariables.AnalyteType;

        /// <summary>
        /// The unmatched-modification tests below pin the process-wide analyte type and write databases of their
        /// own; both are undone here rather than at the end of each test, because a failing assertion would skip
        /// the latter and leave the static set for everything that runs after it.
        /// </summary>
        [TearDown]
        public void RestoreProcessWideState()
        {
            GlobalVariables.AnalyteType = _analyteTypeBeforeTest;
            if (!Directory.Exists(ScratchDirectory))
                return;

            Directory.Delete(ScratchDirectory, true);
            // also take the shared parent, so nothing of this fixture's is left in the test output root
            string scratchRoot = Path.GetDirectoryName(ScratchDirectory);
            if (!Directory.EnumerateFileSystemEntries(scratchRoot).Any())
                Directory.Delete(scratchRoot);
        }

        private const string TestPeptideSequence = "MPSPEPTIDEK";
        private const string TestOligoSequence = "GUCCUGCCUCUAGUGAAGCA";

        /// <summary>
        /// Writes a minimal protein database with one "modified residue" feature per description, all on one entry.
        /// Deliberately omits the &lt;modification&gt; blocks a MetaMorpheus-written database embeds, because those
        /// are loaded by GetPtmListFromProteinXml and would make any description resolvable.
        /// </summary>
        private static string WriteDbWithModificationDescriptions(string fileName, params string[] modificationDescriptions)
        {
            return WriteDb(fileName, TestPeptideSequence, modificationDescriptions);
        }

        /// <summary>
        /// As above, with one &lt;entry&gt; per element of <paramref name="descriptionsPerEntry"/>, so one database can
        /// repeat an unmatched id across entries or carry more unmatched ids than the warning names individually.
        /// Every feature is annotated on position 3 — a serine in the peptide sequence, a cytosine in the oligo one.
        /// </summary>
        private static string WriteDb(string fileName, string sequence, params string[][] descriptionsPerEntry)
        {
            string entries = string.Concat(descriptionsPerEntry.Select((descriptions, entryIndex) =>
                "  <entry>\n" +
                $"    <accession>P1234{entryIndex}</accession>\n" +
                $"    <name>TEST_ENTRY_{entryIndex}</name>\n" +
                "    <protein><recommendedName><fullName>Test entry</fullName></recommendedName></protein>\n" +
                string.Concat(descriptions.Select(description =>
                    $"    <feature type=\"modified residue\" description=\"{description}\">\n" +
                    "      <location>\n" +
                    "        <position position=\"3\" />\n" +
                    "      </location>\n" +
                    "    </feature>\n")) +
                $"    <sequence length=\"{sequence.Length}\">{sequence}</sequence>\n" +
                "  </entry>\n"));

            Directory.CreateDirectory(ScratchDirectory);
            string path = Path.Combine(ScratchDirectory, fileName);
            File.WriteAllText(path,
                "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n" +
                "<mzLibProteinDb>\n" +
                entries +
                "</mzLibProteinDb>\n");
            return path;
        }

        private static List<string> LoadAndCollectWarnings(params string[] dbPaths)
        {
            DatabaseLoadingEngine.LoadBioPolymers("TestTask", dbPaths.Select(p => new DbForTask(p, false)).ToList(),
                true, DecoyType.None, new List<string> { "UniProt" }, new CommonParameters(), out var errors);
            return errors;
        }

        /// <summary>
        /// A modification annotated in a database but absent from the known modifications is skipped by mzLib, which
        /// records it in its unknownModifications dictionary. That dictionary used to be discarded here, so the
        /// annotation vanished with nothing said about it (mzLib #417). It must now reach the user as a warning.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public void UnmatchedModificationInDatabaseIsWarnedAbout()
        {
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
            string dbPath = WriteDbWithModificationDescriptions("unknownModDb.xml", "Nonexistent modification on S");

            var errors = LoadAndCollectWarnings(dbPath);

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
        /// could be satisfied by warning unconditionally, which would be worse than the silence it replaces. The
        /// loaded protein is inspected as well, because an empty unknownModifications dictionary is also what a
        /// database that never parsed produces, and that would satisfy the absence assertion on its own.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public void ResolvableModificationProducesNoUnmatchedModificationWarning()
        {
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
            string dbPath = WriteDbWithModificationDescriptions("knownModDb.xml", "Phosphoserine on S");

            var bioPolymers = DatabaseLoadingEngine.LoadBioPolymers("TestTask",
                new List<DbForTask> { new DbForTask(dbPath, false) },
                true, DecoyType.None, new List<string> { "UniProt" }, new CommonParameters(), out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(errors.Any(p => p.Contains("could not be matched")), Is.False,
                    $"no unmatched-modification warning expected; got: {string.Join(" | ", errors)}");
                Assert.That(bioPolymers, Has.Count.EqualTo(1), "the test database should have yielded one protein");
                Assert.That(bioPolymers[0].OneBasedPossibleLocalizedModifications.TryGetValue(3, out var modsAtThree)
                            && modsAtThree.Any(m => m.IdWithMotif == "Phosphoserine on S"),
                    "the annotation must have been parsed and matched, not merely absent from the unknown set");
            });
        }

        /// <summary>
        /// A database annotating more unmatched modifications than the warning names individually: the first five
        /// in ordinal order are named and the remainder is reported as a count, so a database using many unknown
        /// ids produces one bounded message rather than an unbounded list.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public void UnmatchedModificationWarningNamesFiveThenCountsTheRest()
        {
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
            string[] descriptions = Enumerable.Range(1, 7)
                .Select(i => $"Nonexistent modification {i} on S")
                .ToArray();
            string dbPath = WriteDbWithModificationDescriptions("manyUnknownModsDb.xml", descriptions);

            var errors = LoadAndCollectWarnings(dbPath);

            string warning = errors.SingleOrDefault(p => p.Contains("could not be matched"));
            Assert.That(warning, Is.Not.Null,
                $"expected exactly one unmatched-modification warning; got: {string.Join(" | ", errors)}");
            Assert.Multiple(() =>
            {
                // the count covers every unmatched id, not just the named ones
                Assert.That(warning, Does.Contain("7 distinct annotated modification type(s)"));
                Assert.That(warning, Does.Contain("and 2 more"));
                Assert.That(warning, Does.Contain("'Nonexistent modification 1 on S'"));
                Assert.That(warning, Does.Contain("'Nonexistent modification 5 on S'"));
                Assert.That(warning, Does.Not.Contain("'Nonexistent modification 6 on S'"));
            });
        }

        /// <summary>
        /// Aggregation is per database, which is the reason the warning is emitted inside the per-database loop and
        /// names a file: two databases must produce two warnings, each naming only its own file and its own ids.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public void UnmatchedModificationWarningsAreReportedPerDatabase()
        {
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
            string firstDb = WriteDbWithModificationDescriptions("firstUnknownModDb.xml", "Nonexistent modification A on S");
            string secondDb = WriteDbWithModificationDescriptions("secondUnknownModDb.xml", "Nonexistent modification B on S");

            var warnings = LoadAndCollectWarnings(firstDb, secondDb)
                .Where(p => p.Contains("could not be matched")).ToList();

            Assert.That(warnings, Has.Count.EqualTo(2),
                $"expected one warning per database; got: {string.Join(" | ", warnings)}");
            string firstWarning = warnings.SingleOrDefault(p => p.Contains("firstUnknownModDb.xml"));
            string secondWarning = warnings.SingleOrDefault(p => p.Contains("secondUnknownModDb.xml"));
            Assert.Multiple(() =>
            {
                Assert.That(firstWarning, Is.Not.Null,
                    $"expected a warning naming firstUnknownModDb.xml; got: {string.Join(" | ", warnings)}");
                Assert.That(secondWarning, Is.Not.Null,
                    $"expected a warning naming secondUnknownModDb.xml; got: {string.Join(" | ", warnings)}");
            });
            Assert.Multiple(() =>
            {
                Assert.That(firstWarning, Does.Contain("'Nonexistent modification A on S'"));
                Assert.That(firstWarning, Does.Not.Contain("'Nonexistent modification B on S'"));
                Assert.That(secondWarning, Does.Contain("'Nonexistent modification B on S'"));
                Assert.That(secondWarning, Does.Not.Contain("'Nonexistent modification A on S'"));
            });
        }

        /// <summary>
        /// The other half of per-database aggregation: one unmatched id annotated on three entries is one warning
        /// reporting one modification type, not one warning per entry.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public void RepeatedUnmatchedModificationIsWarnedAboutOnce()
        {
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
            string[] oneUnknownMod = { "Nonexistent modification on S" };
            string dbPath = WriteDb("repeatedUnknownModDb.xml", TestPeptideSequence,
                oneUnknownMod, oneUnknownMod, oneUnknownMod);

            var warnings = LoadAndCollectWarnings(dbPath)
                .Where(p => p.Contains("could not be matched")).ToList();

            Assert.That(warnings, Has.Count.EqualTo(1),
                $"expected a single aggregated warning; got: {string.Join(" | ", warnings)}");
            Assert.That(warnings[0], Does.Contain("1 distinct annotated modification type(s)"));
        }

        /// <summary>
        /// The oligo call site is wired separately from the protein one, so it needs its own coverage: the same
        /// warning must come out of RnaDbLoader.LoadRnaXML.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public void UnmatchedModificationInOligoDatabaseIsWarnedAbout()
        {
            GlobalVariables.AnalyteType = AnalyteType.Oligo;
            string dbPath = WriteDb("unknownOligoModDb.xml", TestOligoSequence,
                new[] { "Nonexistent modification on C" });

            DatabaseLoadingEngine.LoadBioPolymers("TestTask", new List<DbForTask> { new DbForTask(dbPath, false) },
                true, DecoyType.None, new List<string> { "UniProt" },
                new CommonParameters(digestionParams: new RnaDigestionParams()), out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(errors.Any(p => p.Contains("Nonexistent modification on C")),
                    $"expected the unmatched modification to be named; got: {string.Join(" | ", errors)}");
                Assert.That(errors.Any(p => p.Contains("unknownOligoModDb.xml")),
                    "expected the warning to name the database it came from");
            });
        }

        /// <summary>
        /// FASTA is the only input that leaves unknownModifications null, so it is the only thing exercising the
        /// null guard. Dropping that guard would throw on every FASTA load.
        /// </summary>
        [Test]
        [NonParallelizable] // mutates the process-wide GlobalVariables.AnalyteType
        public void FastaDatabaseProducesNoUnmatchedModificationWarning()
        {
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
            string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "gapdh.fasta");

            var bioPolymers = DatabaseLoadingEngine.LoadBioPolymers("TestTask",
                new List<DbForTask> { new DbForTask(dbPath, false) },
                true, DecoyType.None, new List<string> { "UniProt" }, new CommonParameters(), out var errors);

            Assert.Multiple(() =>
            {
                Assert.That(bioPolymers, Is.Not.Empty, "the FASTA should have loaded");
                Assert.That(errors.Any(p => p.Contains("could not be matched")), Is.False,
                    $"no unmatched-modification warning expected from a FASTA; got: {string.Join(" | ", errors)}");
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