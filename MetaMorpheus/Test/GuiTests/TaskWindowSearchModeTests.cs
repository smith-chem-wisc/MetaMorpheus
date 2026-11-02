using GuiFunctions.Util;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using TaskLayer;
using Transcriptomics.Digestion;

namespace Test.GuiTests
{
    /// <summary>
    /// What the task windows load and save for <c>SearchModeType</c> and <c>FragmentationTerminus</c>.
    /// </summary>
    /// <remarks>
    /// <para><b>The bug this guards.</b> The Glyco, crosslink, GPTMD and calibration windows built their DigestionParams
    /// without a search mode or terminus, so saving a task reset both to Full and Both. A semi-specific glyco task (for
    /// example a settings file naming the removed "semi-trypsin", which now loads as trypsin + Semi) silently became fully
    /// specific the moment someone opened and saved it. The windows also showed <c>DigestionParams.Protease</c>, which for a
    /// non-specific search is singleN or singleC rather than the protease the user chose.</para>
    /// <para><b>What the windows do now.</b> The Glyco window has a semi-specific choice (Full or Semi, always with terminus
    /// Both, which asks for peptides rather than seeds). The crosslink, GPTMD and calibration windows have no such control
    /// and keep whatever the loaded task had; the run-time check refuses a combination those searches cannot use. The WPF
    /// windows cannot be referenced from this test project, so the rules live in <see cref="TaskWindowSearchMode"/>.</para>
    /// </remarks>
    [TestFixture]
    public static class TaskWindowSearchModeTests
    {
        [Test]
        [TestCase(CleavageSpecificity.Full, FragmentationTerminus.Both, false)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.Both, true)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.N, true)]
        [TestCase(CleavageSpecificity.None, FragmentationTerminus.N, false)]
        public static void IsSemiSpecific_IsTrueExactlyForASemiSearchMode(CleavageSpecificity searchModeType, FragmentationTerminus terminus, bool expected)
        {
            var loaded = new DigestionParams("trypsin", searchModeType: searchModeType, fragmentationTerminus: terminus);
            Assert.That(TaskWindowSearchMode.IsSemiSpecific(loaded), Is.EqualTo(expected));
        }

        /// <summary>
        /// The Glyco window's choice always saves terminus Both: a glyco search needs the peptides themselves, and N or C
        /// would ask mzLib for seeds only the non-specific search engine can use.
        /// </summary>
        [Test]
        [TestCase(true, CleavageSpecificity.Semi)]
        [TestCase(false, CleavageSpecificity.Full)]
        public static void ForSemiSpecificChoice_SavesThatSearchModeWithBothTermini(bool semiSpecific, CleavageSpecificity expected)
        {
            var (searchModeType, terminus) = TaskWindowSearchMode.ForSemiSpecificChoice(semiSpecific);
            Assert.That(searchModeType, Is.EqualTo(expected));
            Assert.That(terminus, Is.EqualTo(FragmentationTerminus.Both));
        }

        /// <summary>
        /// A glyco settings file naming "semi-trypsin" loads as trypsin + Semi; opening it in the Glyco window and saving
        /// without touching anything keeps it semi-specific.
        /// </summary>
        [Test]
        [NonParallelizable] // loading raises a warning through the static MetaMorpheusTask.WarnHandler
        public static void GlycoWindow_LoadingAndSavingALegacySemiTrypsinTask_KeepsItSemiSpecific()
        {
            var task = MetaMorpheusTask.ReadTaskTomlWithLowResFallback<GlycoSearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "GlycoSearchTaskconfig_ETD.toml"));
            var loaded = (DigestionParams)task.CommonParameters.DigestionParams;

            bool checkBoxShown = TaskWindowSearchMode.IsSemiSpecific(loaded);
            var (searchModeType, terminus) = TaskWindowSearchMode.ForSemiSpecificChoice(checkBoxShown);
            var saved = new DigestionParams(TaskWindowSearchMode.ProteaseToShow(loaded).Name, loaded.MaxMissedCleavages,
                searchModeType: searchModeType, fragmentationTerminus: terminus);

            Assert.That(saved.SearchModeType, Is.EqualTo(CleavageSpecificity.Semi));
            Assert.That(saved.FragmentationTerminus, Is.EqualTo(FragmentationTerminus.Both));
            Assert.That(saved.Protease, Is.SameAs(loaded.Protease));
        }

        /// <summary>
        /// The crosslink, GPTMD and calibration windows keep the loaded search mode and terminus, whatever they are; a new
        /// or RNA task (which has no protein search mode) saves the defaults.
        /// </summary>
        [Test]
        [TestCase(CleavageSpecificity.Full, FragmentationTerminus.Both)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.Both)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.C)]
        [TestCase(CleavageSpecificity.None, FragmentationTerminus.N)]
        public static void Preserve_KeepsTheLoadedSearchModeAndTerminus(CleavageSpecificity searchModeType, FragmentationTerminus terminus)
        {
            var loaded = new DigestionParams("trypsin", searchModeType: searchModeType, fragmentationTerminus: terminus);
            Assert.That(TaskWindowSearchMode.Preserve(loaded), Is.EqualTo((searchModeType, terminus)));
        }

        [Test]
        public static void Preserve_ForATaskWithoutProteinDigestionParams_SavesTheDefaults()
        {
            Assert.That(TaskWindowSearchMode.Preserve(null), Is.EqualTo((CleavageSpecificity.Full, FragmentationTerminus.Both)));
            Assert.That(TaskWindowSearchMode.Preserve(new RnaDigestionParams()), Is.EqualTo((CleavageSpecificity.Full, FragmentationTerminus.Both)));
        }

        /// <summary>
        /// A task the windows can show as it is gets no warning. That covers every task a window itself can save, plus RNA
        /// tasks and tasks without digestion parameters, which have no protein search mode.
        /// </summary>
        [Test]
        [TestCase(CleavageSpecificity.Full, FragmentationTerminus.Both)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.Both)]
        public static void Warnings_ForAUsableSearchMode_AreNull(CleavageSpecificity searchModeType, FragmentationTerminus terminus)
        {
            var loaded = new DigestionParams("trypsin", searchModeType: searchModeType, fragmentationTerminus: terminus);
            Assert.That(TaskWindowSearchMode.ForSemiSpecificChoiceWarning(loaded), Is.Null);
            Assert.That(TaskWindowSearchMode.PreserveWarning(loaded), Is.Null);
        }

        [Test]
        public static void Warnings_ForATaskWithoutProteinDigestionParams_AreNull()
        {
            foreach (IDigestionParams loaded in new IDigestionParams[] { null, new RnaDigestionParams() })
            {
                Assert.That(TaskWindowSearchMode.ForSemiSpecificChoiceWarning(loaded), Is.Null);
                Assert.That(TaskWindowSearchMode.PreserveWarning(loaded), Is.Null);
            }
        }

        private static IEnumerable<TestCaseData> SeedRequests()
        {
            yield return new TestCaseData(CleavageSpecificity.Semi, FragmentationTerminus.N);
            yield return new TestCaseData(CleavageSpecificity.Semi, FragmentationTerminus.C);
            yield return new TestCaseData(CleavageSpecificity.None, FragmentationTerminus.N);
            yield return new TestCaseData(CleavageSpecificity.None, FragmentationTerminus.Both);
        }

        /// <summary>
        /// Review of #2812 (nbollis): a hand-edited glyco task with Semi + N, or None, was rewritten to Semi + Both or
        /// Full + Both on save without a word, so the user searched something other than what the file asked for. Saving
        /// still writes a mode a glyco search can use, but the window now says, before the user saves, which setting was
        /// loaded and that saving replaces it with the choice in the semi-specific box.
        /// </summary>
        [Test]
        [TestCaseSource(nameof(SeedRequests))]
        public static void ForSemiSpecificChoiceWarning_ForASeedRequest_NamesTheLoadedSettingAndSaysSavingChangesIt(CleavageSpecificity searchModeType, FragmentationTerminus terminus)
        {
            var loaded = new DigestionParams("trypsin", searchModeType: searchModeType, fragmentationTerminus: terminus);
            string warning = TaskWindowSearchMode.ForSemiSpecificChoiceWarning(loaded);

            Assert.That(warning, Does.Contain($"SearchModeType {searchModeType}"));
            Assert.That(warning, Does.Contain($"FragmentationTerminus {terminus}"));
            Assert.That(warning, Does.Contain("Saving will change it"));
            Assert.That(warning, Does.Contain("Semi-specific digestion"), "it must point at the box that decides what is saved");
        }

        /// <summary>
        /// Review of #2812 (nbollis): the crosslink, GPTMD and calibration windows keep a loaded mode they have no control
        /// for, so the user saw an ordinary protease and only learned at run time that the task would be refused. The window
        /// now says so up front, and says the fix has to be made in the settings file.
        /// </summary>
        [Test]
        [TestCaseSource(nameof(SeedRequests))]
        public static void PreserveWarning_ForASeedRequest_NamesTheLoadedSettingAndSaysTheTaskWillBeRefused(CleavageSpecificity searchModeType, FragmentationTerminus terminus)
        {
            var loaded = new DigestionParams("trypsin", searchModeType: searchModeType, fragmentationTerminus: terminus);
            string warning = TaskWindowSearchMode.PreserveWarning(loaded);

            Assert.That(warning, Does.Contain($"SearchModeType {searchModeType}"));
            Assert.That(warning, Does.Contain($"FragmentationTerminus {terminus}"));
            Assert.That(warning, Does.Contain("refused"));
            Assert.That(warning, Does.Contain("settings file"));
        }

        /// <summary>
        /// The windows warn about exactly the tasks the run-time check refuses (for a task other than the non-specific
        /// search), so the two can never disagree.
        /// </summary>
        [Test]
        public static void Warnings_AreGivenExactlyWhenTheRunTimeCheckRefusesTheTask()
        {
            foreach (CleavageSpecificity searchModeType in new[] { CleavageSpecificity.Full, CleavageSpecificity.Semi, CleavageSpecificity.None })
            foreach (FragmentationTerminus terminus in new[] { FragmentationTerminus.Both, FragmentationTerminus.N, FragmentationTerminus.C })
            {
                var loaded = new DigestionParams("trypsin", searchModeType: searchModeType, fragmentationTerminus: terminus);
                var task = new GptmdTask { CommonParameters = new EngineLayer.CommonParameters(digestionParams: loaded) };
                bool refused = task.GetSeedDigestionRefusal() != null;

                Assert.That(TaskWindowSearchMode.PreserveWarning(loaded) != null, Is.EqualTo(refused), $"{searchModeType} + {terminus}");
                Assert.That(TaskWindowSearchMode.ForSemiSpecificChoiceWarning(loaded) != null, Is.EqualTo(refused), $"{searchModeType} + {terminus}");
            }
        }

        /// <summary>
        /// For a non-specific search, DigestionParams.Protease is singleN or singleC; the window must show, and save, the
        /// protease the user chose, which is SpecificProtease. Saving singleN as the protease would lose trypsin's cleavage
        /// sites, which still limit missed cleavages.
        /// </summary>
        [Test]
        [TestCase(CleavageSpecificity.Full, FragmentationTerminus.Both)]
        [TestCase(CleavageSpecificity.None, FragmentationTerminus.N)]
        [TestCase(CleavageSpecificity.None, FragmentationTerminus.C)]
        public static void ProteaseToShow_IsTheProteaseTheUserChose(CleavageSpecificity searchModeType, FragmentationTerminus terminus)
        {
            var loaded = new DigestionParams("trypsin", searchModeType: searchModeType, fragmentationTerminus: terminus);
            Assert.That(TaskWindowSearchMode.ProteaseToShow(loaded).Name, Is.EqualTo("trypsin"));

            var resaved = new DigestionParams(TaskWindowSearchMode.ProteaseToShow(loaded).Name,
                searchModeType: TaskWindowSearchMode.Preserve(loaded).SearchModeType, fragmentationTerminus: TaskWindowSearchMode.Preserve(loaded).Terminus);
            Assert.That(resaved, Is.EqualTo(loaded), "loading and saving without changes must give back the same digestion parameters");
        }
    }
}
