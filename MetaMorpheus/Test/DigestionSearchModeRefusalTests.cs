using EngineLayer;
using EngineLayer.DatabaseLoading;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using Transcriptomics.Digestion;

namespace Test
{
    /// <summary>
    /// A task whose digestion settings ask for seed peptides, but which cannot use seeds, is refused with a message
    /// instead of silently searching the wrong peptides.
    /// </summary>
    /// <remarks>
    /// <para><b>Peptides versus seeds.</b> In mzLib, <c>SearchModeType</c> Full gives fully specific peptides and Semi with
    /// <c>FragmentationTerminus</c> Both gives semi-specific peptides. Semi with N or C, and None with any terminus, give
    /// seeds: long stretches fixed at one terminus whose other end is decided after the search, from the precursor mass.
    /// Only the non-specific search engine (a Search task with SearchType NonSpecific) does that. Every other search
    /// (Classic, Modern, Glyco, crosslink, GPTMD, calibration) scores what digestion gives it as it is, so with seeds it
    /// runs to completion and quietly reports far fewer, wrong identifications. That is how semi-specific glyco searches
    /// lost most of their identifications before mzLib #1303.</para>
    /// <para><b>Where the refusal happens.</b> <see cref="EverythingRunnerEngine"/> checks every task before running any of
    /// them, so hours of earlier tasks are not wasted, and refuses with a warning, following the non-specific-RNA refusal
    /// (#2759): an exception out of a task reaches the GUI as a crash, and the user never sees the message saying what to
    /// change. <see cref="MetaMorpheusTask.RunTask"/> throws the same message as a backstop for callers that bypass the
    /// runner.</para>
    /// </remarks>
    [TestFixture]
    public static class DigestionSearchModeRefusalTests
    {
        private static CommonParameters Common(CleavageSpecificity searchModeType, FragmentationTerminus terminus) =>
            new(digestionParams: new DigestionParams("trypsin", searchModeType: searchModeType, fragmentationTerminus: terminus));

        private static SearchTask Search(SearchType searchType, CleavageSpecificity searchModeType, FragmentationTerminus terminus) =>
            new() { CommonParameters = Common(searchModeType, terminus), SearchParameters = new SearchParameters { SearchType = searchType } };

        /// <summary>
        /// Builds the task inside the test, not in the case source: some task constructors (crosslink) need global settings
        /// that are only loaded once the test run has started.
        /// </summary>
        private static MetaMorpheusTask MakeTask(string kind, CleavageSpecificity searchModeType, FragmentationTerminus terminus) => kind switch
        {
            "Glyco" => new GlycoSearchTask { CommonParameters = Common(searchModeType, terminus) },
            "Crosslink" => new XLSearchTask { CommonParameters = Common(searchModeType, terminus) },
            "Gptmd" => new GptmdTask { CommonParameters = Common(searchModeType, terminus) },
            "Calibration" => new CalibrationTask { CommonParameters = Common(searchModeType, terminus) },
            "ClassicSearch" => Search(SearchType.Classic, searchModeType, terminus),
            "ModernSearch" => Search(SearchType.Modern, searchModeType, terminus),
            "NonSpecificSearch" => Search(SearchType.NonSpecific, searchModeType, terminus),
            "GlycoWithRnaDigestion" => new GlycoSearchTask { CommonParameters = new CommonParameters(digestionParams: new RnaDigestionParams()) },
            _ => throw new System.ArgumentException(kind),
        };

        private static IEnumerable<TestCaseData> EveryTaskAndDigestion()
        {
            foreach (string kind in new[] { "Glyco", "Crosslink", "Gptmd", "Calibration", "ClassicSearch", "ModernSearch" })
            {
                yield return new TestCaseData(kind, CleavageSpecificity.Full, FragmentationTerminus.Both, false).SetName($"{kind} Full+Both is allowed");
                yield return new TestCaseData(kind, CleavageSpecificity.Semi, FragmentationTerminus.Both, false).SetName($"{kind} Semi+Both is allowed (semi-specific peptides)");
                yield return new TestCaseData(kind, CleavageSpecificity.Semi, FragmentationTerminus.N, true).SetName($"{kind} Semi+N is refused (seeds)");
                yield return new TestCaseData(kind, CleavageSpecificity.Semi, FragmentationTerminus.C, true).SetName($"{kind} Semi+C is refused (seeds)");
                yield return new TestCaseData(kind, CleavageSpecificity.None, FragmentationTerminus.N, true).SetName($"{kind} None+N is refused (seeds)");
                yield return new TestCaseData(kind, CleavageSpecificity.None, FragmentationTerminus.Both, true).SetName($"{kind} None+Both is refused (seeds)");
            }

            // the non-specific search trims seeds, so every combination the Search window can save is allowed
            foreach (var (mode, terminus) in new[] { (CleavageSpecificity.Semi, FragmentationTerminus.N), (CleavageSpecificity.Semi, FragmentationTerminus.C),
                         (CleavageSpecificity.Semi, FragmentationTerminus.Both), (CleavageSpecificity.None, FragmentationTerminus.N),
                         (CleavageSpecificity.None, FragmentationTerminus.C), (CleavageSpecificity.None, FragmentationTerminus.Both),
                         (CleavageSpecificity.Full, FragmentationTerminus.Both) })
                yield return new TestCaseData("NonSpecificSearch", mode, terminus, false).SetName($"NonSpecificSearch {mode}+{terminus} is allowed");

            // RNA digestion has no protein seed request
            yield return new TestCaseData("GlycoWithRnaDigestion", CleavageSpecificity.Full, FragmentationTerminus.Both, false).SetName("RNA digestion parameters are not checked");
        }

        /// <summary>The rule for every task type and every combination of search mode and terminus.</summary>
        [Test]
        [TestCaseSource(nameof(EveryTaskAndDigestion))]
        public static void GetRefusal_RefusesExactlyTheSeedRequestsOfTasksThatCannotUseSeeds(string taskKind, CleavageSpecificity searchModeType, FragmentationTerminus terminus, bool expectedRefused)
        {
            string refusal = DigestionSearchModeCheck.GetRefusal(MakeTask(taskKind, searchModeType, terminus));
            if (expectedRefused)
            {
                Assert.That(refusal, Does.StartWith("Cannot proceed."));
                Assert.That(refusal, Does.Contain("seed"), "the message must say why: these settings give seeds");
                Assert.That(refusal, Does.Contain("FragmentationTerminus Both").Or.Contain("non-specific search"),
                    "the message must say what to do instead");
            }
            else
            {
                Assert.That(refusal, Is.Null);
            }
        }

        /// <summary>
        /// The runner checks every task before it runs any: a valid first task does not run when a later task is refused,
        /// the refusal reaches the user as a warning, and nothing is thrown (a throw would reach the GUI as a crash).
        /// </summary>
        [Test]
        [NonParallelizable] // EverythingRunnerEngine warnings go through static handlers
        public static void Runner_WithOneTaskThatCannotUseItsSeeds_RunsNoTaskAndWarns()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, nameof(Runner_WithOneTaskThatCannotUseItsSeeds_RunsNoTaskAndWarns));
            if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            Directory.CreateDirectory(outputFolder);
            try
            {
                var tasks = new List<(string, MetaMorpheusTask)>
                {
                    ("ValidClassicSearch", Search(SearchType.Classic, CleavageSpecificity.Full, FragmentationTerminus.Both)),
                    ("GlycoWithSeeds", new GlycoSearchTask { CommonParameters = Common(CleavageSpecificity.Semi, FragmentationTerminus.N) }),
                };
                var runner = new EverythingRunnerEngine(tasks,
                    new List<string> { Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf") },
                    new List<DbForTask> { new(Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "P16150.fasta"), false) },
                    outputFolder);

                Assert.DoesNotThrow(() => runner.Run(), "the runner refuses; it does not fault a task");
                Assert.That(runner.Warnings.Count(w => w.StartsWith("Cannot proceed.") && w.Contains("GlycoWithSeeds")), Is.EqualTo(1),
                    "the refusal names the task and reaches the user as a warning");
                Assert.That(Directory.Exists(Path.Combine(outputFolder, "ValidClassicSearch")), Is.False, "no task runs, not even the valid one before it");
                Assert.That(Directory.Exists(Path.Combine(outputFolder, "GlycoWithSeeds")), Is.False);
            }
            finally
            {
                if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>A caller that runs a task directly, bypassing the runner, gets the same message as an exception.</summary>
        [Test]
        public static void RunTask_CalledDirectlyWithSeedsItCannotUse_ThrowsTheRefusal()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, nameof(RunTask_CalledDirectlyWithSeedsItCannotUse_ThrowsTheRefusal));
            string taskFolder = Path.Combine(outputFolder, "Task");
            if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            Directory.CreateDirectory(taskFolder);
            Directory.CreateDirectory(Path.Combine(outputFolder, "Task Settings"));
            try
            {
                var task = new GlycoSearchTask { CommonParameters = Common(CleavageSpecificity.None, FragmentationTerminus.N) };
                var thrown = Assert.Throws<MetaMorpheusException>(() => task.RunTask(taskFolder,
                    new List<DbForTask> { new(Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "P16150.fasta"), false) },
                    new List<string> { Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf") },
                    "Task"));
                Assert.That(thrown.Message, Is.EqualTo(DigestionSearchModeCheck.GetRefusal(task, "Task")));
            }
            finally
            {
                if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            }
        }
    }
}
