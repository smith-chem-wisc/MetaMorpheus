using EngineLayer;
using MetaMorpheusCommandLine;
using Nett;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using Transcriptomics;
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
            "Averaging" => new SpectralAveragingTask { CommonParameters = Common(searchModeType, terminus) },
            _ => throw new System.ArgumentException(kind),
        };

        /// <summary>A task of this kind that digests nucleic acids with this rnase, search mode and terminus.</summary>
        private static MetaMorpheusTask MakeRnaTask(string kind, string rnase,
            CleavageSpecificity searchModeType = CleavageSpecificity.Full, FragmentationTerminus terminus = FragmentationTerminus.Both)
        {
            var common = new CommonParameters(digestionParams: new RnaDigestionParams(rnase, fragmentationTerminus: terminus, searchModeType: searchModeType));
            return kind switch
            {
                "Gptmd" => new GptmdTask { CommonParameters = common },
                "Calibration" => new CalibrationTask { CommonParameters = common },
                "ClassicSearch" => new SearchTask { CommonParameters = common, SearchParameters = new SearchParameters { SearchType = SearchType.Classic } },
                "ModernSearch" => new SearchTask { CommonParameters = common, SearchParameters = new SearchParameters { SearchType = SearchType.Modern } },
                "NonSpecificSearch" => new SearchTask { CommonParameters = common, SearchParameters = new SearchParameters { SearchType = SearchType.NonSpecific } },
                "Averaging" => new SpectralAveragingTask { CommonParameters = common },
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

            // a spectral averaging task digests nothing, so no digestion settings are refused for it
            foreach (var (mode, terminus) in new[] { (CleavageSpecificity.Semi, FragmentationTerminus.N), (CleavageSpecificity.Semi, FragmentationTerminus.C),
                         (CleavageSpecificity.None, FragmentationTerminus.Both), (CleavageSpecificity.Full, FragmentationTerminus.Both) })
                yield return new TestCaseData("Averaging", mode, terminus, false).SetName($"Averaging {mode}+{terminus} is allowed (it does not digest)");
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
        /// With SearchModeType Full the rnase decides: singleN and singleC give seed oligos, every other rnase gives oligos.
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

        /// <summary>
        /// Since mzLib 1.0.592 nucleic acid digestion follows SearchModeType: None with any terminus, and Semi with a single
        /// terminus, give seed oligos from a fully specific rnase. Semi with Both asks for semi-specific oligos, not seeds.
        /// </summary>
        [Test]
        [TestCase(CleavageSpecificity.None, FragmentationTerminus.Both, true)]
        [TestCase(CleavageSpecificity.None, FragmentationTerminus.FivePrime, true)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.FivePrime, true)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.ThreePrime, true)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.N, true)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.C, true)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.Both, false)]
        [TestCase(CleavageSpecificity.Full, FragmentationTerminus.FivePrime, false)]
        public static void AsksForSeeds_ForNucleicAcids_FollowsSearchModeType(CleavageSpecificity searchModeType, FragmentationTerminus terminus, bool expected)
        {
            var digestionParams = new RnaDigestionParams("RNase T1", fragmentationTerminus: terminus, searchModeType: searchModeType);
            Assert.That(MetaMorpheusTask.AsksForSeeds(digestionParams), Is.EqualTo(expected));
        }

        /// <summary>A task without CommonParameters has no digestion settings to refuse.</summary>
        [Test]
        public static void GetSeedDigestionRefusal_WithoutCommonParameters_IsNull()
        {
            var task = new SearchTask { CommonParameters = null };
            Assert.That(task.GetSeedDigestionRefusal("MyTask"), Is.Null);
        }

        private static IEnumerable<TestCaseData> EveryRnaTaskAndRnase()
        {
            const CleavageSpecificity Full = CleavageSpecificity.Full, Semi = CleavageSpecificity.Semi, None = CleavageSpecificity.None;
            const FragmentationTerminus Both = FragmentationTerminus.Both;
            foreach (string kind in new[] { "Gptmd", "Calibration", "ClassicSearch", "ModernSearch" })
            {
                yield return new TestCaseData(kind, "RNase T1", Full, Both, null).SetName($"{kind} with RNase T1 is allowed");
                yield return new TestCaseData(kind, "top-down", Full, Both, null).SetName($"{kind} with top-down is allowed");
                yield return new TestCaseData(kind, "non-specific", Full, Both, null).SetName($"{kind} with the rnase non-specific is allowed");
                yield return new TestCaseData(kind, "singleN", Full, Both, "the rnase singleN").SetName($"{kind} with singleN is refused (seed oligos)");
                yield return new TestCaseData(kind, "singleC", Full, Both, "the rnase singleC").SetName($"{kind} with singleC is refused (seed oligos)");
                yield return new TestCaseData(kind, "RNase T1", None, Both, "SearchModeType None").SetName($"{kind} with RNase T1 and SearchModeType None is refused (seed oligos)");
                yield return new TestCaseData(kind, "RNase T1", Semi, FragmentationTerminus.ThreePrime, "SearchModeType Semi with FragmentationTerminus ThreePrime")
                    .SetName($"{kind} with RNase T1, Semi and ThreePrime is refused (seed oligos)");
            }

            // non-specific search over nucleic acids is refused by the runner with its own message, not by this check
            yield return new TestCaseData("NonSpecificSearch", "singleN", Full, Both, null).SetName("NonSpecificSearch with singleN is left to the runner");
            yield return new TestCaseData("NonSpecificSearch", "RNase T1", Semi, FragmentationTerminus.FivePrime, null).SetName("NonSpecificSearch with Semi and FivePrime is left to the runner");

            // a spectral averaging task digests nothing, so its rnase is never used
            yield return new TestCaseData("Averaging", "singleC", Full, Both, null).SetName("Averaging with singleC is allowed (it does not digest)");
        }

        /// <summary>
        /// A nucleic acid task that asks for seed oligos is refused, and the message names the setting and says what to choose
        /// instead. Before this, such a task ran to completion and scored each seed as if it were a whole oligo.
        /// </summary>
        [Test]
        [TestCaseSource(nameof(EveryRnaTaskAndRnase))]
        public static void GetSeedDigestionRefusal_ForNucleicAcids_RefusesSeedDigestion(string taskKind, string rnase,
            CleavageSpecificity searchModeType, FragmentationTerminus terminus, string expectedSetting)
        {
            string refusal = MakeRnaTask(taskKind, rnase, searchModeType, terminus).GetSeedDigestionRefusal("MyTask");
            if (expectedSetting != null)
            {
                Assert.That(refusal, Does.StartWith("Cannot proceed."));
                Assert.That(refusal, Does.Contain("\"MyTask\""), "the message must name the task");
                Assert.That(refusal, Does.Contain(expectedSetting), "the message must name the setting that asks for seeds");
                Assert.That(refusal, Does.Contain("seed oligos"), "the message must say why");
                Assert.That(refusal, Does.Contain("fully specific rnase"), "the message must say what to do instead");
            }
            else
            {
                Assert.That(refusal, Is.Null);
            }
        }

        /// <summary>
        /// mzLib takes Semi from the rnase's own specificity as well as from SearchModeType, so an rnase defined as semi (a
        /// user can add one in rnase_custom.tsv) gives seed oligos with a single terminus even under SearchModeType Full. Such
        /// a task is refused like RNase T1 with SearchModeType Semi, and the message names the rnase. With Both it is left
        /// alone, as SearchModeType Semi with Both is.
        /// </summary>
        [Test]
        [TestCase(FragmentationTerminus.FivePrime)]
        [TestCase(FragmentationTerminus.ThreePrime)]
        public static void GetSeedDigestionRefusal_WithASemiSpecificRnaseAndOneTerminus_RefusesTheTask(FragmentationTerminus terminus)
        {
            const string semiRnase = "semi G| (DigestionSearchModeRefusalTests)";
            string tsv = Path.Combine(TestContext.CurrentContext.WorkDirectory, $"{nameof(GetSeedDigestionRefusal_WithASemiSpecificRnaseAndOneTerminus_RefusesTheTask)}_{terminus}.tsv");
            File.WriteAllLines(tsv, new[] { "Name\tMotif\tSpecificity", $"{semiRnase}\tG|\tsemi" });
            try
            {
                Assert.That(RnaseDictionary.LoadAndMergeCustomRnases(tsv).Added, Is.EqualTo(new[] { semiRnase }));

                // premise: under Full, mzLib digests the semi rnase into the seeds of RNase T1 with SearchModeType Semi
                var rna = new RNA("AAGUACUGCCGAUUGAGGCAU");
                var semiRnaseUnderFull = new RnaDigestionParams(semiRnase, maxMissedCleavages: 1, fragmentationTerminus: terminus);
                var rnaseT1UnderSemi = new RnaDigestionParams("RNase T1", maxMissedCleavages: 1, fragmentationTerminus: terminus, searchModeType: CleavageSpecificity.Semi);
                Assert.That(Oligos(rna, semiRnaseUnderFull), Is.EquivalentTo(Oligos(rna, rnaseT1UnderSemi)));
                Assert.That(MetaMorpheusTask.AsksForSeeds(rnaseT1UnderSemi), Is.True);

                Assert.That(MetaMorpheusTask.AsksForSeeds(semiRnaseUnderFull), Is.True);
                foreach (string kind in new[] { "Gptmd", "Calibration", "ClassicSearch", "ModernSearch" })
                {
                    string refusal = MakeRnaTask(kind, semiRnase, CleavageSpecificity.Full, terminus).GetSeedDigestionRefusal("MyTask");
                    Assert.That(refusal, Does.Contain($"uses the semi-specific rnase {semiRnase} with FragmentationTerminus {terminus}"), kind);
                    Assert.That(refusal, Does.Contain("seed oligos"), kind);
                }

                Assert.That(MetaMorpheusTask.AsksForSeeds(new RnaDigestionParams(semiRnase)), Is.False, "with Both it is semi-specific, not seeds");
            }
            finally
            {
                RnaseDictionary.Dictionary.Remove(semiRnase);
                File.Delete(tsv);
            }
        }

        private static List<string> Oligos(RNA rna, RnaDigestionParams digestionParams) =>
            rna.Digest(digestionParams, new List<Modification>(), new List<Modification>()).Select(o => o.BaseSequence).Distinct().ToList();

        /// <summary>
        /// A task file saved before mzLib 1.0.592 has no SpecificRnase, so that reads back as the parameterless constructor's
        /// top-down while Rnase says singleN or singleC. The refusal is decided on Rnase, so the message has to name Rnase,
        /// not fall through to a Semi wording for a file that says SearchModeType Full.
        /// </summary>
        [Test]
        [TestCase("singleN")]
        [TestCase("singleC")]
        public static void GetSeedDigestionRefusal_ForATaskFileWithoutSpecificRnase_NamesTheRnase(string rnase)
        {
            string text = File.ReadAllText(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "RnaSearchTask_MultipleDecon.toml"));
            Assert.That(text, Does.Not.Contain("SpecificRnase"), "premise: the file predates SpecificRnase");
            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, $"{nameof(GetSeedDigestionRefusal_ForATaskFileWithoutSpecificRnase_NamesTheRnase)}_{rnase}.toml");
            File.WriteAllText(path, text.Replace("Rnase = \"top-down\"", $"Rnase = \"{rnase}\""));
            try
            {
                var task = Toml.ReadFile<SearchTask>(path, MetaMorpheusTask.tomlConfig);
                var digestionParams = (RnaDigestionParams)task.CommonParameters.DigestionParams;
                Assert.That(digestionParams.Rnase.Name, Is.EqualTo(rnase), "premise");
                Assert.That(digestionParams.SearchModeType, Is.EqualTo(CleavageSpecificity.Full), "premise");

                string refusal = task.GetSeedDigestionRefusal("Task1SearchTask");
                Assert.That(refusal, Does.Contain($"uses the rnase {rnase},"));
                Assert.That(refusal, Does.Not.Contain("SearchModeType Semi"));
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// Writes a file-specific toml next to a spectra file in its own folder and returns the spectra file's path. The
        /// spectra file itself need not exist: only the toml beside it is read before a run.
        /// </summary>
        private static string SpectraFileWithToml(string testName, string toml)
        {
            string folder = Path.Combine(TestContext.CurrentContext.WorkDirectory, testName);
            Directory.CreateDirectory(folder);
            string spectraFile = Path.Combine(folder, "sample.mzML");
            File.WriteAllText(Path.Combine(folder, "sample.toml"), toml);
            return spectraFile;
        }

        /// <summary>
        /// Each spectra file is digested with what <see cref="MetaMorpheusTask.SetAllFileSpecificCommonParams"/> builds from
        /// its file-specific toml, and for nucleic acids that takes the rnase from the file. A file that names singleN or
        /// singleC is therefore searched as seed oligos even when the task's own rnase is fully specific, so the refusal must
        /// check every file's digestion, and name the file.
        /// </summary>
        [Test]
        [TestCase("Gptmd", "DigestionAgent = \"singleC\"", "uses the rnase singleC")]
        [TestCase("Calibration", "DigestionAgent = \"singleN\"", "uses the rnase singleN")]
        [TestCase("ClassicSearch", "Rnase = \"singleC\"", "uses the rnase singleC")]
        [TestCase("ModernSearch", "DigestionAgent = \"singleN\"", "uses the rnase singleN")]
        [TestCase("ClassicSearch", "DigestionAgent = \"RNase T1\"", null)]
        [TestCase("ClassicSearch", "DigestionAgent = \"no such rnase\"", null)] // unreadable: the run warns about it, as before
        [TestCase("NonSpecificSearch", "DigestionAgent = \"singleC\"", null)]
        public static void GetSeedDigestionRefusal_WithAFileSpecificRnaseThatAsksForSeeds_RefusesAndNamesTheFile(string taskKind, string toml, string expectedSetting)
        {
            string spectraFile = SpectraFileWithToml($"{nameof(GetSeedDigestionRefusal_WithAFileSpecificRnaseThatAsksForSeeds_RefusesAndNamesTheFile)}_{TestContext.CurrentContext.Test.ID}", toml);
            MetaMorpheusTask task = MakeRnaTask(taskKind, "RNase T1");
            Assert.That(task.GetSeedDigestionRefusal("MyTask"), Is.Null, "premise: the task's own rnase is fully specific");

            string refusal = task.GetSeedDigestionRefusal("MyTask", new List<string> { spectraFile });
            if (expectedSetting == null)
            {
                Assert.That(refusal, Is.Null);
                return;
            }
            Assert.That(refusal, Does.StartWith("Cannot proceed."));
            Assert.That(refusal, Does.Contain("\"MyTask\""));
            Assert.That(refusal, Does.Contain("file-specific parameters of sample.mzML"), "the message must name the file");
            Assert.That(refusal, Does.Contain(expectedSetting));
            Assert.That(refusal, Does.Contain("seed oligos"));
        }

        /// <summary>
        /// A custom semi rnase named in a file-specific toml gives seed oligos when the task's terminus is FivePrime or
        /// ThreePrime, as it does when the task names it.
        /// </summary>
        [Test]
        public static void GetSeedDigestionRefusal_WithAFileSpecificSemiRnaseAndOneTerminus_Refuses()
        {
            const string semiRnase = "semi G| per file (DigestionSearchModeRefusalTests)";
            string tsv = Path.Combine(TestContext.CurrentContext.WorkDirectory, $"{nameof(GetSeedDigestionRefusal_WithAFileSpecificSemiRnaseAndOneTerminus_Refuses)}.tsv");
            File.WriteAllLines(tsv, new[] { "Name\tMotif\tSpecificity", $"{semiRnase}\tG|\tsemi" });
            try
            {
                Assert.That(RnaseDictionary.LoadAndMergeCustomRnases(tsv).Added, Is.EqualTo(new[] { semiRnase }));
                string spectraFile = SpectraFileWithToml(nameof(GetSeedDigestionRefusal_WithAFileSpecificSemiRnaseAndOneTerminus_Refuses), $"DigestionAgent = \"{semiRnase}\"");
                MetaMorpheusTask task = MakeRnaTask("ClassicSearch", "RNase T1", CleavageSpecificity.Full, FragmentationTerminus.FivePrime);
                Assert.That(task.GetSeedDigestionRefusal("MyTask"), Is.Null, "premise");

                string refusal = task.GetSeedDigestionRefusal("MyTask", new List<string> { spectraFile });
                Assert.That(refusal, Does.Contain($"uses the semi-specific rnase {semiRnase} with FragmentationTerminus FivePrime"));
            }
            finally
            {
                RnaseDictionary.Dictionary.Remove(semiRnase);
                File.Delete(tsv);
            }
        }

        [Test]
        public static void CommandLine_WithAFileSpecificRnaseThatAsksForSeeds_RefusesTheRun()
        {
            string spectraFile = SpectraFileWithToml(nameof(CommandLine_WithAFileSpecificRnaseThatAsksForSeeds_RefusesTheRun), "DigestionAgent = \"singleC\"");
            var tasks = new List<(string, MetaMorpheusTask)> { ("Task1SearchTask", MakeRnaTask("ClassicSearch", "RNase T1")) };
            var output = new StringWriter();

            Assert.That(Program.RefuseSeedDigestion(tasks, output, true), Is.EqualTo(0), "premise: without the spectra files it carries on");
            Assert.That(Program.RefuseSeedDigestion(tasks, output, true, new List<string> { spectraFile }), Is.EqualTo(6));
            Assert.That(output.ToString().Trim(), Is.EqualTo(tasks[0].Item2.GetSeedDigestionRefusal("Task1SearchTask", new List<string> { spectraFile })));
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
