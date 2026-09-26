using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Truncation;
using MassSpectrometry;
using MzLibUtil;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
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
            // Locked defaults (docs/Truncation-Search.md #3, #17).
            Assert.That(task.TruncationSearchParameters.ParentQValueThreshold, Is.EqualTo(0.10));
            Assert.That(task.TruncationSearchParameters.WriteDecoys, Is.True);
            Assert.That(task.TruncationSearchParameters.WriteContaminants, Is.True);
        }

        /// <summary>Only acceptors that map a monoisotopic precursor to a chop at fixed notches are allowed (#9).</summary>
        [TestCase(MassDiffAcceptorType.Exact, true)]
        [TestCase(MassDiffAcceptorType.ThreeMM, true)]
        [TestCase(MassDiffAcceptorType.PlusOrMinusThreeMM, true)]
        [TestCase(MassDiffAcceptorType.Open, false)]
        [TestCase(MassDiffAcceptorType.ModOpen, false)]
        [TestCase(MassDiffAcceptorType.MostAbundant_Exact, false)]
        [TestCase(MassDiffAcceptorType.MostAbundant_PlusMinusOne, false)]
        public void ValidateChopAcceptor_RejectsAcceptorsThatCannotChop(MassDiffAcceptorType type, bool allowed)
        {
            MassDiffAcceptor acceptor = SearchTask.GetMassDiffAcceptor(new PpmTolerance(10), type, null);
            if (allowed)
                Assert.DoesNotThrow(() => TruncationSearchTask.ValidateChopAcceptor(acceptor, type));
            else
                Assert.Throws<MetaMorpheusException>(() => TruncationSearchTask.ValidateChopAcceptor(acceptor, type));
        }

        [Test]
        public void DescribeParameterDifferences_ReportsOnlyTheSettingsThatDiffer()
        {
            var search = new CommonParameters(precursorMassTolerance: new PpmTolerance(10), productMassTolerance: new PpmTolerance(20));
            var same = new CommonParameters(precursorMassTolerance: new PpmTolerance(10), productMassTolerance: new PpmTolerance(20));
            var looser = new CommonParameters(precursorMassTolerance: new PpmTolerance(5), productMassTolerance: new PpmTolerance(20));

            Assert.That(TruncationSearchTask.DescribeParameterDifferences(search, same), Is.Empty);
            List<string> differences = TruncationSearchTask.DescribeParameterDifferences(search, looser);
            Assert.That(differences, Has.Count.EqualTo(1));
            Assert.That(differences[0], Does.StartWith("PrecursorMassTolerance"));
        }

        /// <summary>WriteDecoys / WriteContaminants decide which rows reach the output files (#17).</summary>
        [TestCase(true, true, 3)]
        [TestCase(false, true, 2)]
        [TestCase(true, false, 2)]
        [TestCase(false, false, 1)]
        public void RowsToWrite_HonorsWriteDecoysAndWriteContaminants(bool writeDecoys, bool writeContaminants, int expectedRows)
        {
            var cp = new CommonParameters();
            var scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new[] { 100.0 }, new[] { 1.0 }, false), 1, 2, true,
                Polarity.Positive, 1, new MzRange(0, 1000), "f", MZAnalyzerType.Orbitrap, 1, 1, null, "scan=1"), 500, 1, "f", cp);
            SpectralMatch Psm(string accession, bool isDecoy, bool isContaminant)
            {
                var protein = new Protein("PEPTIDE", accession, isDecoy: isDecoy, isContaminant: isContaminant);
                var form = protein.Digest(new DigestionParams(protease: "top-down"), new List<Modification>(), new List<Modification>()).First();
                var psm = new PeptideSpectralMatch(form, 0, 10, 0, scan, cp, new List<Omics.Fragmentation.MatchedFragmentIon>());
                psm.ResolveAllAmbiguities();
                return psm;
            }

            var parameters = new TruncationSearchParameters { WriteDecoys = writeDecoys, WriteContaminants = writeContaminants };
            List<SpectralMatch> rows = TruncationSearchTask.RowsToWrite(
                new[] { Psm("T", false, false), Psm("DECOY_T", true, false), Psm("C", false, true), null }, parameters);

            Assert.That(rows, Has.Count.EqualTo(expectedRows));
            Assert.That(rows.Any(p => p.IsDecoy), Is.EqualTo(writeDecoys));
            Assert.That(rows.Any(p => p.IsContaminant), Is.EqualTo(writeContaminants));
        }

        /// <summary>
        /// A disk parent reported at [2 to 120] (Met cleaved) keeps those coordinates, so a 5-residue C-chop is
        /// (2-115), not (1-114), and a parent that starts with M away from residue 1 is not read as NME (#13).
        /// </summary>
        [Test]
        public void DiskParent_KeepsProteinCoordinatesThroughTheChop()
        {
            string sequence = "ADEKRHSTNQGVLIFYWPC" + "ADEKRHSTNQGVLIFYWPC";
            var parent = TruncationSearchTask.BuildDiskProteoform(sequence, "P", 2, false, new DigestionParams());
            Assert.That((parent.OneBasedStartResidueInProtein, parent.OneBasedEndResidueInProtein), Is.EqualTo((2, 39)));
            var exact = SearchTask.GetMassDiffAcceptor(new PpmTolerance(10), MassDiffAcceptorType.Exact, null);

            double cChopMass = TruncationSearchTask.BuildDiskProteoform(sequence.Substring(0, sequence.Length - 5), "x", 1, false, new DigestionParams()).MonoisotopicMass;
            ChopResult cChop = ProteoformChopper.ChopUntilMassMatches(parent, Omics.Fragmentation.FragmentationTerminus.C, cChopMass, exact);
            Assert.That(cChop.TruncatedForm.Description, Is.EqualTo(TruncationPass3.CTerminalTruncation + "(2-34)"));

            var startsWithMet = TruncationSearchTask.BuildDiskProteoform("M" + sequence, "Q", 30, false, new DigestionParams());
            ChopResult nChop = ProteoformChopper.ChopUntilMassMatches(startsWithMet, Omics.Fragmentation.FragmentationTerminus.N, parent.MonoisotopicMass, exact);
            Assert.That(nChop.TruncatedForm.Description, Is.EqualTo(TruncationPass3.NTerminalTruncation + "(31-68)"));
        }

        [Test]
        public void ParseStartResidues_ReadsEachAlternative()
        {
            Assert.That(TruncationSearchTask.ParseStartResidues("[2 to 120]"), Is.EqualTo(new[] { 2 }));
            Assert.That(TruncationSearchTask.ParseStartResidues("[2 to 120]|[31 to 149]"), Is.EqualTo(new[] { 2, 31 }));
            Assert.That(TruncationSearchTask.ParseStartResidues(null), Is.Empty);
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

        /// <summary>
        /// The perf log is benchmarking instrumentation, not a search setting: it must not appear in the
        /// TOML every user sees, and a TOML that still carries the key from before it moved must still
        /// load. The environment variable is the supported way to turn it on.
        /// </summary>
        [Test]
        public void PerfLogPath_IsNotPartOfTheTomlSurface()
        {
            var task = new TruncationSearchTask();
            task.TruncationSearchParameters.PerfLogPath = @"C:\somewhere\perf_log.tsv";

            string tomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TruncationNoPerfLog.toml");
            Toml.WriteFile(task, tomlPath, MetaMorpheusTask.tomlConfig);
            try
            {
                Assert.That(File.ReadAllText(tomlPath), Does.Not.Contain("PerfLogPath"));
                Assert.That(File.ReadAllText(tomlPath), Does.Not.Contain("CustomOutputFolderName"));

                // A stale key left over in an existing benchmarking TOML must not break the read.
                File.AppendAllText(tomlPath, "\r\nPerfLogPath = 'C:\\stale\\perf_log.tsv'\r\n");
                var read = Toml.ReadFile<TruncationSearchTask>(tomlPath, MetaMorpheusTask.tomlConfig);
                Assert.That(read.TruncationSearchParameters.PerfLogPath, Is.Null,
                    "the TOML key is ignored; the environment variable drives it");
            }
            finally
            {
                File.Delete(tomlPath);
            }
        }

        /// <summary>The benchmarking hook is driven by an environment variable, not by the task TOML.</summary>
        [Test]
        public void PerfLogPath_DefaultsFromEnvironmentVariable()
        {
            string previous = Environment.GetEnvironmentVariable(TruncationSearchParameters.PerfLogPathEnvironmentVariable);
            try
            {
                Environment.SetEnvironmentVariable(TruncationSearchParameters.PerfLogPathEnvironmentVariable, @"C:\bench\perf_log.tsv");
                Assert.That(new TruncationSearchParameters().PerfLogPath, Is.EqualTo(@"C:\bench\perf_log.tsv"));

                Environment.SetEnvironmentVariable(TruncationSearchParameters.PerfLogPathEnvironmentVariable, null);
                Assert.That(new TruncationSearchParameters().PerfLogPath, Is.Null);
            }
            finally
            {
                Environment.SetEnvironmentVariable(TruncationSearchParameters.PerfLogPathEnvironmentVariable, previous);
            }
        }

        /// <summary>
        /// The runner must not hand out — or hold on to — the task-chain context for run lists that have no
        /// consumer, and must drop the deposited results once a consumer has run.
        /// </summary>
        [Test]
        public void TaskChainContext_OnlyWiredWhenARunListHasAConsumer()
        {
            Assert.That(new SearchTask().ConsumesTaskChainContext, Is.False);
            Assert.That(new TruncationSearchTask().ConsumesTaskChainContext, Is.True);

            var context = new TaskChainContext();
            context.Deposit("Task1", new List<SpectralMatch>());
            Assert.That(context.TryGetMostRecent(out List<SpectralMatch> _), Is.True);

            context.Clear();
            Assert.That(context.TryGetMostRecent(out List<SpectralMatch> _), Is.False);
            Assert.That(context.TryGet("Task1", out List<SpectralMatch> _), Is.False);
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
                // Also exercise the optional perf-log path (docs/Truncation-Search.md, Benchmarking hook).
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

                // Every confident Pass 1 PSM is inherited exactly once as full-length, on its own precursor, with
                // the same Essential Sequence as AllPSMs; none of those precursors is also reported as a truncation (#4a).
                var pass1 = ReadTsv(Directory.GetFiles(Path.Combine(outDirectory, "Task1-SearchTask"), "AllPSMs.psmtsv", SearchOption.AllDirectories).Single());
                var truncation = ReadTsv(psmsPath);
                var confident = pass1.Where(r => r["Decoy/Contaminant/Target"] == "T" && double.Parse(r["QValue"]) <= 0.01).ToList();
                Assert.That(confident, Is.Not.Empty);
                foreach (var row in confident)
                {
                    var samePrecursor = truncation.Where(t => t["Scan Number"] == row["Scan Number"]
                        && Math.Abs(double.Parse(t["Precursor Mass"]) - double.Parse(row["Precursor Mass"])) < 0.01).ToList();
                    Assert.That(samePrecursor.Count(t => t["Description"] == TruncationPass3.FullLength), Is.EqualTo(1),
                        $"scan {row["Scan Number"]} should be inherited once as full-length");
                    Assert.That(samePrecursor, Has.Count.EqualTo(1), $"scan {row["Scan Number"]} also carries a truncation");
                    Assert.That(samePrecursor[0]["Essential Sequence"], Is.EqualTo(row["Essential Sequence"]));
                }

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
        /// proteoforms). Also exercises the perf-log + CandidateRanks output.
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
                task.TruncationSearchParameters.PerfLogPath = Path.Combine(outDir, "perf_log.tsv"); // -> perf log + CandidateRanks

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

        /// <summary>Rows of a psmtsv as column-name -> cell dictionaries.</summary>
        private static List<Dictionary<string, string>> ReadTsv(string path)
        {
            string[] lines = File.ReadAllLines(path);
            string[] header = lines[0].Split('\t');
            return lines.Skip(1)
                .Select(line => line.Split('\t'))
                .Select(cells => header.Select((name, i) => (name, cell: i < cells.Length ? cells[i] : ""))
                    .GroupBy(c => c.name).ToDictionary(g => g.Key, g => g.First().cell))
                .ToList();
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
