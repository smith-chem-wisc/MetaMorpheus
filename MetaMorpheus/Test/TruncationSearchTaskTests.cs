using System.Collections.Generic;
using System.IO;
using EngineLayer;
using EngineLayer.DatabaseLoading;
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
        /// Phase 0.7 gate: a run list of [SearchTask, TruncationSearchTask] executes end-to-end through
        /// EverythingRunnerEngine without throwing, and the (stub) truncation task produces its own
        /// output folder with no truncation TSVs yet. Uses the existing tiny top-down fixture.
        /// </summary>
        [Test]
        public void EverythingRunner_SearchThenTruncation_RunsWithoutThrowing()
        {
            string outDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TruncationScaffoldE2E");
            if (Directory.Exists(outDirectory))
                Directory.Delete(outDirectory, true);

            string msAlignPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "JurkatTopDownRep2Fract1_ms2.msalign");
            string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "ThreeHumanHistone.fasta");
            string topDownSearchToml = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownTestData", "TopDownSearchToml.toml");

            var searchTask = Toml.ReadFile<SearchTask>(topDownSearchToml, MetaMorpheusTask.tomlConfig);
            searchTask.CommonParameters.PrecursorMassTolerance = new AbsoluteTolerance(5);
            var truncationTask = new TruncationSearchTask();

            var taskList = new List<(string, MetaMorpheusTask)>
            {
                ("Task1-SearchTask", searchTask),
                ("Task2-TruncationSearchTask", truncationTask)
            };

            var engine = new EverythingRunnerEngine(taskList, new List<string> { msAlignPath },
                new List<DbForTask> { new DbForTask(dbPath, false) }, outDirectory);

            Assert.DoesNotThrow(() => engine.Run());

            string truncationOutput = Path.Combine(outDirectory, "Task2-TruncationSearchTask");
            Assert.That(Directory.Exists(truncationOutput), Is.True, "TruncationSearchTask did not produce an output folder.");
            Assert.That(Directory.GetFiles(truncationOutput, "AllTruncated*.psmtsv"), Is.Empty, "Stub should emit no truncation TSVs yet.");

            Directory.Delete(outDirectory, true);
        }
    }
}
