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
        public void Scaffolding_TestProjectWired()
        {
            Assert.That(true, Is.True);
        }

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

            Directory.Delete(outDirectory, true);
        }
    }
}
