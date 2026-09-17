using EngineLayer;
using MetaMorpheusCommandLine;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
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
    /// <para><b>Where the refusal happens.</b> Where every other run-level check does, before the run starts: the GUI's
    /// Run button (TaskValidator.CheckDigestionSearchMode, called from MainWindow.RunAllTasks_Click) and the command line
    /// (<see cref="Program.RefuseSeedDigestion"/>, next to the experimental-design check). Tasks themselves assume their
    /// settings are valid. The rule and the message live on <see cref="MetaMorpheusTask"/> so both front ends say the same
    /// thing; the WPF check cannot be referenced from this test project, so the rule and the command line are tested here.</para>
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
            _ => throw new System.ArgumentException(kind),
        };

        /// <summary>A task of this kind that digests nucleic acids with this rnase.</summary>
        private static MetaMorpheusTask MakeRnaTask(string kind, string rnase)
        {
            var common = new CommonParameters(digestionParams: new RnaDigestionParams(rnase));
            return kind switch
            {
                "Gptmd" => new GptmdTask { CommonParameters = common },
                "Calibration" => new CalibrationTask { CommonParameters = common },
                "ClassicSearch" => new SearchTask { CommonParameters = common, SearchParameters = new SearchParameters { SearchType = SearchType.Classic } },
                "ModernSearch" => new SearchTask { CommonParameters = common, SearchParameters = new SearchParameters { SearchType = SearchType.Modern } },
                "NonSpecificSearch" => new SearchTask { CommonParameters = common, SearchParameters = new SearchParameters { SearchType = SearchType.NonSpecific } },
                _ => throw new System.ArgumentException(kind),
            };
        }

        /// <summary>
        /// The one rule for which digestion settings give seeds. The refusal and the task windows' warnings both use it.
        /// </summary>
        [Test]
        [TestCase(CleavageSpecificity.Full, FragmentationTerminus.Both, false)]
        [TestCase(CleavageSpecificity.Full, FragmentationTerminus.N, false)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.Both, false)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.N, true)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.C, true)]
        [TestCase(CleavageSpecificity.None, FragmentationTerminus.Both, true)]
        [TestCase(CleavageSpecificity.None, FragmentationTerminus.C, true)]
        public static void AsksForSeeds_IsTrueForNoneAndForSemiWithOneTerminus(CleavageSpecificity searchModeType, FragmentationTerminus terminus, bool expected)
        {
            var digestionParams = new DigestionParams("trypsin", searchModeType: searchModeType, fragmentationTerminus: terminus);
            Assert.That(MetaMorpheusTask.AsksForSeeds(digestionParams), Is.EqualTo(expected));
        }

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

        }

        /// <summary>The rule for every task type and every combination of search mode and terminus.</summary>
        [Test]
        [TestCaseSource(nameof(EveryTaskAndDigestion))]
        public static void GetSeedDigestionRefusal_RefusesExactlyTheSeedRequestsOfTasksThatCannotUseSeeds(string taskKind, CleavageSpecificity searchModeType, FragmentationTerminus terminus, bool expectedRefused)
        {
            string refusal = MakeTask(taskKind, searchModeType, terminus).GetSeedDigestionRefusal("MyTask");
            if (expectedRefused)
            {
                Assert.That(refusal, Does.StartWith("Cannot proceed."));
                Assert.That(refusal, Does.Contain("\"MyTask\""), "the message must name the task");
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
        /// Nucleic acid digestion ignores SearchModeType; the rnase decides. singleN and singleC give seed oligos, every other
        /// rnase gives oligos.
        /// </summary>
        [Test]
        [TestCase("singleN", true)]
        [TestCase("singleC", true)]
        [TestCase("RNase T1", false)]
        [TestCase("non-specific", false)]
        [TestCase("top-down", false)]
        public static void AsksForSeeds_ForNucleicAcids_IsTrueForSingleNAndSingleC(string rnase, bool expected)
        {
            Assert.That(MetaMorpheusTask.AsksForSeeds(new RnaDigestionParams(rnase)), Is.EqualTo(expected));
        }

        /// <summary>SearchModeType is not what makes a nucleic acid task ask for seeds, so it does not cause a refusal.</summary>
        [Test]
        public static void AsksForSeeds_ForNucleicAcids_IgnoresSearchModeType()
        {
            var digestionParams = new RnaDigestionParams("RNase T1") { SearchModeType = CleavageSpecificity.None };
            Assert.That(MetaMorpheusTask.AsksForSeeds(digestionParams), Is.False);
        }

        private static IEnumerable<TestCaseData> EveryRnaTaskAndRnase()
        {
            foreach (string kind in new[] { "Gptmd", "Calibration", "ClassicSearch", "ModernSearch" })
            {
                yield return new TestCaseData(kind, "RNase T1", false).SetName($"{kind} with RNase T1 is allowed");
                yield return new TestCaseData(kind, "top-down", false).SetName($"{kind} with top-down is allowed");
                yield return new TestCaseData(kind, "non-specific", false).SetName($"{kind} with the rnase non-specific is allowed");
                yield return new TestCaseData(kind, "singleN", true).SetName($"{kind} with singleN is refused (seed oligos)");
                yield return new TestCaseData(kind, "singleC", true).SetName($"{kind} with singleC is refused (seed oligos)");
            }

            // the non-specific search trims seed oligos, so it is the one search that needs them
            yield return new TestCaseData("NonSpecificSearch", "singleN", false).SetName("NonSpecificSearch with singleN is allowed");
            yield return new TestCaseData("NonSpecificSearch", "singleC", false).SetName("NonSpecificSearch with singleC is allowed");
        }

        /// <summary>
        /// A nucleic acid task that asks for seed oligos is refused, and the message says what to choose instead. Before
        /// this, such a task ran to completion and scored each seed as if it were a whole oligo.
        /// </summary>
        [Test]
        [TestCaseSource(nameof(EveryRnaTaskAndRnase))]
        public static void GetSeedDigestionRefusal_ForNucleicAcids_RefusesSingleNAndSingleC(string taskKind, string rnase, bool expectedRefused)
        {
            string refusal = MakeRnaTask(taskKind, rnase).GetSeedDigestionRefusal("MyTask");
            if (expectedRefused)
            {
                Assert.That(refusal, Does.StartWith("Cannot proceed."));
                Assert.That(refusal, Does.Contain("\"MyTask\""), "the message must name the task");
                Assert.That(refusal, Does.Contain(rnase), "the message must name the rnase");
                Assert.That(refusal, Does.Contain("seed oligos"), "the message must say why");
                Assert.That(refusal, Does.Contain("non-specific search type").And.Contain("fully specific rnase"), "the message must say what to do instead");
            }
            else
            {
                Assert.That(refusal, Is.Null);
            }
        }

        /// <summary>
        /// The other direction: the non-specific search over nucleic acids can only trim seeds, so any other rnase is
        /// refused, and so are complementary ions, which have no nucleic acid definition.
        /// </summary>
        [Test]
        [TestCase("RNase T1")]
        [TestCase("top-down")]
        [TestCase("non-specific")]
        public static void GetSeedDigestionRefusal_NonSpecificSearchOverNucleicAcids_NeedsSingleNOrSingleC(string rnase)
        {
            string refusal = MakeRnaTask("NonSpecificSearch", rnase).GetSeedDigestionRefusal("MyTask");

            Assert.That(refusal, Does.StartWith("Cannot proceed."));
            Assert.That(refusal, Does.Contain("\"MyTask\"").And.Contain(rnase));
            Assert.That(refusal, Does.Contain("singleN or singleC").And.Contain("Classic or Modern"));
        }

        [Test]
        public static void GetSeedDigestionRefusal_NonSpecificSearchOverNucleicAcids_RefusesComplementaryIons()
        {
            var task = new SearchTask
            {
                CommonParameters = new CommonParameters(addCompIons: true, digestionParams: new RnaDigestionParams("singleN")),
                SearchParameters = new SearchParameters { SearchType = SearchType.NonSpecific },
            };

            string refusal = task.GetSeedDigestionRefusal("MyTask");

            Assert.That(refusal, Does.StartWith("Cannot proceed."));
            Assert.That(refusal, Does.Contain("complementary ions"));
        }

        [Test]
        public static void CommandLine_WithANucleicAcidTaskThatAsksForSeeds_RefusesTheRun()
        {
            var tasks = new List<(string, MetaMorpheusTask)> { ("Task1SearchTask", MakeRnaTask("ClassicSearch", "singleC")) };
            var output = new StringWriter();

            Assert.That(Program.RefuseSeedDigestion(tasks, output, true), Is.EqualTo(6));
            Assert.That(output.ToString().Trim(), Is.EqualTo(tasks[0].Item2.GetSeedDigestionRefusal("Task1SearchTask")));
        }

        /// <summary>
        /// The command line checks every task before the runner starts, as it does for the experimental design: a valid
        /// first task does not run when a later task is refused, the message names the refused task, and the exit code
        /// tells a script that the run was refused.
        /// </summary>
        [Test]
        public static void CommandLine_WithOneTaskThatCannotUseItsSeeds_RefusesTheRunWithAnExitCode()
        {
            var tasks = new List<(string, MetaMorpheusTask)>
            {
                ("Task1SearchTask", Search(SearchType.Classic, CleavageSpecificity.Full, FragmentationTerminus.Both)),
                ("Task2GlycoSearchTask", MakeTask("Glyco", CleavageSpecificity.Semi, FragmentationTerminus.N)),
            };
            var output = new StringWriter();

            int exitCode = Program.RefuseSeedDigestion(tasks, output, true);

            Assert.That(exitCode, Is.EqualTo(6));
            Assert.That(output.ToString().Trim(), Is.EqualTo(tasks[1].Item2.GetSeedDigestionRefusal("Task2GlycoSearchTask")));
        }

        [Test]
        public static void CommandLine_WithOnlyUsableDigestion_CarriesOnSilently()
        {
            var tasks = new List<(string, MetaMorpheusTask)>
            {
                ("Task1SearchTask", Search(SearchType.NonSpecific, CleavageSpecificity.None, FragmentationTerminus.N)),
                ("Task2GlycoSearchTask", MakeTask("Glyco", CleavageSpecificity.Semi, FragmentationTerminus.Both)),
            };
            var output = new StringWriter();

            Assert.That(Program.RefuseSeedDigestion(tasks, output, true), Is.EqualTo(0));
            Assert.That(output.ToString(), Is.Empty);
        }

        /// <summary>With verbosity none the message is not printed, but the run is still refused.</summary>
        [Test]
        public static void CommandLine_Quiet_StillRefuses()
        {
            var tasks = new List<(string, MetaMorpheusTask)> { ("Task1GptmdTask", MakeTask("Gptmd", CleavageSpecificity.None, FragmentationTerminus.C)) };
            var output = new StringWriter();

            Assert.That(Program.RefuseSeedDigestion(tasks, output, false), Is.EqualTo(6));
            Assert.That(output.ToString(), Is.Empty);
        }
    }
}
