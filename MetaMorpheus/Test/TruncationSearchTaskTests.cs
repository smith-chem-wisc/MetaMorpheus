using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Truncation;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Phase 0 scaffolding tests for the TruncationSearchTask: confirms the task and its in-memory
    /// task-chain hand-off compile and behave as wired. Algorithm-level tests arrive in Phases 1-3.
    /// </summary>
    [TestFixture]
    public class TruncationSearchTaskTests
    {
        [Test]
        public void TruncationSearchTask_Constructs_WithDefaultParameters()
        {
            var task = new TruncationSearchTask();

            Assert.That(task.TaskType, Is.EqualTo(MyTask.Truncation));
            Assert.That(task.CommonParameters, Is.Not.Null);
            Assert.That(task.TruncationSearchParameters, Is.Not.Null);
            // Locked defaults (01_Architecture.md #3, #17).
            Assert.That(task.TruncationSearchParameters.ParentQValueThreshold, Is.EqualTo(0.10));
            Assert.That(task.TruncationSearchParameters.WriteDecoys, Is.True);
            Assert.That(task.TruncationSearchParameters.WriteContaminants, Is.True);
        }

        [Test]
        public void TaskChainContext_RoundTripsTypedResult()
        {
            var context = new TaskChainContext();
            var deposited = new List<string> { "proteoformA", "proteoformB" };

            context.Deposit("Task1-SearchTask", deposited);

            Assert.That(context.TryGet<List<string>>("Task1-SearchTask", out var retrieved), Is.True);
            Assert.That(retrieved, Is.SameAs(deposited));
        }

        [Test]
        public void TaskChainContext_Miss_ReturnsFalseAndDefault()
        {
            var context = new TaskChainContext();

            Assert.That(context.TryGet<List<string>>("nonexistent", out var retrieved), Is.False);
            Assert.That(retrieved, Is.Null);
        }

        [Test]
        public void TaskChainContext_WrongType_ReturnsFalse()
        {
            var context = new TaskChainContext();
            context.Deposit("Task1-SearchTask", 42);

            Assert.That(context.TryGet<List<string>>("Task1-SearchTask", out _), Is.False);
        }

        [Test]
        public void TaskChainContext_TryGetMostRecent_ReturnsLatestAssignable()
        {
            var context = new TaskChainContext();
            var first = new List<string> { "first" };
            var second = new List<string> { "second" };

            context.Deposit("Task1-SearchTask", first);
            context.Deposit("Task2-OtherTask", 7);          // wrong type, skipped
            context.Deposit("Task3-SearchTask", second);

            Assert.That(context.TryGetMostRecent<List<string>>(out var latest), Is.True);
            Assert.That(latest, Is.SameAs(second));

            var empty = new TaskChainContext();
            Assert.That(empty.TryGetMostRecent<List<string>>(out var none), Is.False);
            Assert.That(none, Is.Null);
        }

        /// <summary>
        /// The task must round-trip through TOML so CMD can read a run-list TOML (TaskType = "Truncation")
        /// and dispatch it. Confirms the task type and the truncation-specific settings survive.
        /// </summary>
        [Test]
        public void TruncationSearchTask_TomlRoundTrips()
        {
            var task = new TruncationSearchTask();
            task.TruncationSearchParameters.UpstreamSearchTaskId = "Task1SearchTask";
            task.TruncationSearchParameters.ParentQValueThreshold = 0.05;
            task.TruncationSearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.TwoMM;

            string tomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TruncationRoundTrip.toml");
            Toml.WriteFile(task, tomlPath, MetaMorpheusTask.tomlConfig);
            try
            {
                // CMD dispatches on this raw string (Program.cs switch).
                Assert.That(Toml.ReadFile(tomlPath, MetaMorpheusTask.tomlConfig).Get<string>("TaskType"), Is.EqualTo("Truncation"));

                var read = Toml.ReadFile<TruncationSearchTask>(tomlPath, MetaMorpheusTask.tomlConfig);
                Assert.That(read.TaskType, Is.EqualTo(MyTask.Truncation));
                Assert.That(read.TruncationSearchParameters.UpstreamSearchTaskId, Is.EqualTo("Task1SearchTask"));
                Assert.That(read.TruncationSearchParameters.ParentQValueThreshold, Is.EqualTo(0.05));
                Assert.That(read.TruncationSearchParameters.MassDiffAcceptorType, Is.EqualTo(MassDiffAcceptorType.TwoMM));
            }
            finally
            {
                File.Delete(tomlPath);
            }
        }

        [Test]
        public void PerfLogger_ParsesRunFolderConvention()
        {
            var (phase, dataset, label) = PerfLogger.ParseRunFolderName("2026-05-22_Phase4_KaulichSPE_snippetFromFile");
            Assert.That(phase, Is.EqualTo("Phase4"));
            Assert.That(dataset, Is.EqualTo("KaulichSPE"));
            Assert.That(label, Is.EqualTo("snippetFromFile"));

            // Non-conforming names degrade gracefully.
            var (phase2, dataset2, label2) = PerfLogger.ParseRunFolderName("TruncationE2E");
            Assert.That(phase2, Is.EqualTo("Manual"));
            Assert.That(dataset2, Is.EqualTo(""));
            Assert.That(label2, Is.EqualTo(""));
        }

        [Test]
        public void PerfLogger_WritesHeaderOnceThenAppends()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "perf_log_unit.tsv");
            if (File.Exists(path)) File.Delete(path);
            try
            {
                PerfLogger.Append(path, new TruncationPerfMetrics { NPsmsEmitted = 5, NTruncationsNterm = 3 });
                PerfLogger.Append(path, new TruncationPerfMetrics { NPsmsEmitted = 9, NTruncationsCterm = 2 });

                string[] lines = File.ReadAllLines(path);
                Assert.That(lines.Length, Is.EqualTo(3), "Expected one header + two data rows.");
                Assert.That(lines[0].Split('\t'), Does.Contain("n_psms_emitted"));
                int col = System.Array.IndexOf(lines[0].Split('\t'), "n_psms_emitted");
                Assert.That(lines[1].Split('\t')[col], Is.EqualTo("5"));
                Assert.That(lines[2].Split('\t')[col], Is.EqualTo("9"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void PerfLogger_EmptyOrNullFolder_ReturnsManualDefaults()
        {
            foreach (string input in new[] { "", null })
            {
                var (phase, dataset, label) = PerfLogger.ParseRunFolderName(input);
                Assert.That(phase, Is.EqualTo("Manual"));
                Assert.That(dataset, Is.EqualTo(""));
                Assert.That(label, Is.EqualTo(""));
            }
        }

        [Test]
        public void PerfLogger_SanitizesTabsAndNewlinesInStringCells()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "perf_log_sanitize.tsv");
            if (File.Exists(path)) File.Delete(path);
            try
            {
                // A raw tab/newline in a string field would otherwise split into extra columns / extra rows.
                PerfLogger.Append(path, new TruncationPerfMetrics { RunLabel = "a\tb\nc", DatasetTag = "x\ty" });

                string[] lines = File.ReadAllLines(path);
                int headerCols = lines[0].Split('\t').Length;
                Assert.That(lines.Length, Is.EqualTo(2));                          // header + one row (no extra line from \n)
                Assert.That(lines[1].Split('\t').Length, Is.EqualTo(headerCols));  // no extra columns from \t
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// Phase 3.3 gate: a run list of [SearchTask, TruncationSearchTask] executes end-to-end through
        /// EverythingRunnerEngine without throwing, and the truncation task — fed the upstream search's
        /// proteoforms via the in-memory <see cref="TaskChainContext"/> (#1) — writes both result files
        /// with the standard psmtsv header. Uses the existing tiny top-down fixture.
        /// </summary>
        [Test]
        public void EverythingRunner_SearchThenTruncation_WritesWellFormedOutputs()
        {
            string outDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TruncationE2E");
            if (Directory.Exists(outDirectory))
                Directory.Delete(outDirectory, true);

            try
            {
                // Canonical tiny top-down fixture that actually produces proteoform hits (mirrors the
                // EverythingRunnerEngine TopDownQValue case): sliced yeast TD data + small yeast DB.
                string dataPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "slicedTDYeast.mzML");
                string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "smalldb.fasta");
                string topDownSearchToml = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TopDownSearchToml.toml");

                var searchTask = Toml.ReadFile<SearchTask>(topDownSearchToml, MetaMorpheusTask.tomlConfig);

                var truncationTask = new TruncationSearchTask();
                // Consume the upstream search's proteoforms in-memory (#1) and use the same top-down-tuned
                // CommonParameters so MS2 deconvolution + the intact-skip tolerance match Pass 1.
                truncationTask.TruncationSearchParameters.UpstreamSearchTaskId = "Task1-SearchTask";
                truncationTask.CommonParameters = Toml.ReadFile<SearchTask>(topDownSearchToml, MetaMorpheusTask.tomlConfig).CommonParameters;
                // Also exercise the optional perf-log path (03_Benchmarks).
                string perfLogPath = Path.Combine(outDirectory, "perf_log.tsv");
                truncationTask.TruncationSearchParameters.PerfLogPath = perfLogPath;

                var taskList = new List<(string, MetaMorpheusTask)>
                {
                    ("Task1-SearchTask", searchTask),
                    ("Task2-TruncationSearchTask", truncationTask)
                };

                var engine = new EverythingRunnerEngine(taskList, new List<string> { dataPath },
                    new List<DbForTask> { new DbForTask(dbPath, false) }, outDirectory);

                Assert.DoesNotThrow(() => engine.Run());

                string truncationOutput = Path.Combine(outDirectory, "Task2-TruncationSearchTask");
                Assert.That(Directory.Exists(truncationOutput), Is.True, "TruncationSearchTask did not produce an output folder.");

                string psmsPath = Path.Combine(truncationOutput, TruncationSearchTask.TruncatedPsmsFileName);
                string proteoformsPath = Path.Combine(truncationOutput, TruncationSearchTask.TruncatedProteoformsFileName);
                Assert.That(File.Exists(psmsPath), Is.True, "AllTruncatedPSMs.psmtsv was not written.");
                Assert.That(File.Exists(proteoformsPath), Is.True, "AllTruncatedProteoforms.psmtsv was not written.");

                // Both files carry the standard psmtsv header and the pooled run produced at least one row
                // (the upstream search's intact match, inherited as a full-length form per #4a).
                string expectedHeader = SpectralMatch.GetTabSeparatedHeader();
                string[] psmLines = File.ReadAllLines(psmsPath);
                string[] proteoformLines = File.ReadAllLines(proteoformsPath);
                Assert.That(psmLines[0].TrimEnd(), Is.EqualTo(expectedHeader.TrimEnd()));
                Assert.That(proteoformLines[0].TrimEnd(), Is.EqualTo(expectedHeader.TrimEnd()));
                Assert.That(psmLines.Length, Is.GreaterThan(1), "Expected at least one pooled PSM row.");
                Assert.That(proteoformLines.Length, Is.GreaterThan(1), "Expected at least one pooled proteoform row.");

                // Truncation type rides the Description column (#13); the inherited intact match is "full-length" (#4a).
                Assert.That(psmLines.Skip(1).Any(line => line.Contains(TruncationPass3.FullLength)),
                    Is.True, "Expected an inherited full-length row in the truncation PSM output.");

                // Perf log got one TruncationSearchTask row (header + >=1 data row).
                Assert.That(File.Exists(perfLogPath), Is.True, "perf_log.tsv was not written.");
                string[] perfLines = File.ReadAllLines(perfLogPath);
                Assert.That(perfLines.Length, Is.GreaterThanOrEqualTo(2));
                string[] perfHeader = perfLines[0].Split('\t');
                int taskTypeCol = System.Array.IndexOf(perfHeader, "task_type");
                int nPsmsCol = System.Array.IndexOf(perfHeader, "n_psms_emitted");
                string[] perfRow = perfLines[1].Split('\t');
                Assert.That(perfRow[taskTypeCol], Is.EqualTo("TruncationSearchTask"));
                Assert.That(int.Parse(perfRow[nPsmsCol]), Is.EqualTo(psmLines.Length - 1), "perf n_psms_emitted should match the PSM TSV row count.");
            }
            finally
            {
                if (Directory.Exists(outDirectory))
                    Directory.Delete(outDirectory, true);
            }
        }

        // ---------- alternate parent-seeding paths (run the task standalone so there is no in-memory
        //            TaskChainContext hand-off, which forces the database / tag-filtered / disk paths) ----------

        /// <summary>
        /// Database-seeded path: parents are generated by digesting the protein DB (not from upstream
        /// proteoforms). Also exercises the perf-log + CandidateRanks output, the internal-fragment search,
        /// and the UsePerScanTagRestriction-without-prerequisites warning (decision #12).
        /// </summary>
        [Test]
        public void TruncationTask_DatabaseSeeded_WritesWellFormedOutputs()
        {
            var (data, db, toml) = TopDownFixture();
            string outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TruncationDbSeeded");
            if (Directory.Exists(outDir)) Directory.Delete(outDir, true);
            try
            {
                var task = new TruncationSearchTask
                {
                    CommonParameters = Toml.ReadFile<SearchTask>(toml, MetaMorpheusTask.tomlConfig).CommonParameters
                };
                task.TruncationSearchParameters.SeedParentsFromDatabase = true;    // -> BuildParentsFromDatabase
                task.TruncationSearchParameters.MaxParentMass = 100000;            // don't exclude large yeast proteins
                task.TruncationSearchParameters.SearchInternalTruncations = true;  // -> internal-search branch + file
                task.TruncationSearchParameters.PerfLogPath = Path.Combine(outDir, "perf_log.tsv"); // -> perf log + CandidateRanks
                task.TruncationSearchParameters.UsePerScanTagRestriction = true;   // prereqs off -> exercises the #12 warning

                RunTruncationTaskAlone(task, data, db, outDir);
                AssertTruncationOutputsWellFormed(Path.Combine(outDir, "Task1-TruncationSearchTask"));
            }
            finally
            {
                if (Directory.Exists(outDir)) Directory.Delete(outDir, true);
            }
        }

        /// <summary>
        /// Database-seeded + sequence-tag-filtered path: narrow the DB to tag-supported proteins, with the
        /// per-scan restriction on. Exercises BuildParentsFromDatabaseTagFiltered (tag extraction, k-mer
        /// index, candidate union, parallel digest) and the per-scan allowed-accession map.
        /// </summary>
        [Test]
        public void TruncationTask_TagFilteredSeeding_WritesWellFormedOutputs()
        {
            var (data, db, toml) = TopDownFixture();
            string outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TruncationTagFiltered");
            if (Directory.Exists(outDir)) Directory.Delete(outDir, true);
            try
            {
                var task = new TruncationSearchTask
                {
                    CommonParameters = Toml.ReadFile<SearchTask>(toml, MetaMorpheusTask.tomlConfig).CommonParameters
                };
                task.TruncationSearchParameters.SeedParentsFromDatabase = true;
                task.TruncationSearchParameters.UseSequenceTagFilter = true;       // -> BuildParentsFromDatabaseTagFiltered
                task.TruncationSearchParameters.UsePerScanTagRestriction = true;   // -> per-scan allowed map
                task.TruncationSearchParameters.MinTagHits = 1;                    // permissive so candidates are found
                task.TruncationSearchParameters.MaxParentMass = 100000;

                RunTruncationTaskAlone(task, data, db, outDir);
                AssertTruncationOutputsWellFormed(Path.Combine(outDir, "Task1-TruncationSearchTask"));
            }
            finally
            {
                if (Directory.Exists(outDir)) Directory.Delete(outDir, true);
            }
        }

        /// <summary>
        /// Disk-ingest path: first run the search alone to produce an AllProteoforms.psmtsv (a separate
        /// engine, so there is no in-memory hand-off), then point the truncation task at that file to force
        /// BuildParentsFromDisk + the disk parent filter.
        /// </summary>
        [Test]
        public void TruncationTask_DiskSeeded_WritesWellFormedOutputs()
        {
            var (data, db, toml) = TopDownFixture();
            string outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TruncationDiskSeeded");
            if (Directory.Exists(outDir)) Directory.Delete(outDir, true);
            try
            {
                string searchOut = Path.Combine(outDir, "search");
                var searchTask = Toml.ReadFile<SearchTask>(toml, MetaMorpheusTask.tomlConfig);
                new EverythingRunnerEngine(
                    new List<(string, MetaMorpheusTask)> { ("Task1-SearchTask", searchTask) },
                    new List<string> { data }, new List<DbForTask> { new DbForTask(db, false) }, searchOut).Run();

                string proteoformsTsv = Directory
                    .GetFiles(searchOut, "AllProteoforms.psmtsv", SearchOption.AllDirectories)
                    .FirstOrDefault();
                Assert.That(proteoformsTsv, Is.Not.Null, "search did not produce AllProteoforms.psmtsv to seed the disk path");

                string truncOut = Path.Combine(outDir, "trunc");
                var task = new TruncationSearchTask
                {
                    CommonParameters = Toml.ReadFile<SearchTask>(toml, MetaMorpheusTask.tomlConfig).CommonParameters
                };
                task.TruncationSearchParameters.Pass1ProteoformsFilePath = proteoformsTsv; // -> BuildParentsFromDisk

                RunTruncationTaskAlone(task, data, db, truncOut);
                AssertTruncationOutputsWellFormed(Path.Combine(truncOut, "Task1-TruncationSearchTask"));
            }
            finally
            {
                if (Directory.Exists(outDir)) Directory.Delete(outDir, true);
            }
        }

        // ---------- helpers for the standalone-task integration tests ----------

        private static (string data, string db, string toml) TopDownFixture()
        {
            string dir = TestContext.CurrentContext.TestDirectory;
            return (Path.Combine(dir, "TopDownTestData", "slicedTDYeast.mzML"),
                    Path.Combine(dir, "TestData", "smalldb.fasta"),
                    Path.Combine(dir, "TopDownTestData", "TopDownSearchToml.toml"));
        }

        private static void RunTruncationTaskAlone(TruncationSearchTask task, string data, string db, string outDir)
        {
            var taskList = new List<(string, MetaMorpheusTask)> { ("Task1-TruncationSearchTask", task) };
            var engine = new EverythingRunnerEngine(taskList, new List<string> { data },
                new List<DbForTask> { new DbForTask(db, false) }, outDir);
            Assert.DoesNotThrow(() => engine.Run());
        }

        private static void AssertTruncationOutputsWellFormed(string truncationTaskOutDir)
        {
            Assert.That(Directory.Exists(truncationTaskOutDir), Is.True, "truncation task produced no output folder");
            string psms = Path.Combine(truncationTaskOutDir, TruncationSearchTask.TruncatedPsmsFileName);
            string proteoforms = Path.Combine(truncationTaskOutDir, TruncationSearchTask.TruncatedProteoformsFileName);
            Assert.That(File.Exists(psms), Is.True, "AllTruncatedPSMs.psmtsv not written");
            Assert.That(File.Exists(proteoforms), Is.True, "AllTruncatedProteoforms.psmtsv not written");
            string header = SpectralMatch.GetTabSeparatedHeader().TrimEnd();
            Assert.That(File.ReadAllLines(psms)[0].TrimEnd(), Is.EqualTo(header));
            Assert.That(File.ReadAllLines(proteoforms)[0].TrimEnd(), Is.EqualTo(header));
        }
    }
}
