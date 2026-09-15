using GuiFunctions.Util;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
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
