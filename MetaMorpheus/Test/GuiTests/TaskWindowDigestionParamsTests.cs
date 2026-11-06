using GuiFunctions.Util;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Transcriptomics.Digestion;

namespace Test.GuiTests
{
    /// <summary>
    /// What the task windows save for the digestion settings they do not show.
    /// </summary>
    /// <remarks>
    /// <para><b>The bug this guards.</b> Each task window builds a new <see cref="DigestionParams"/> on save from its own
    /// controls. Every constructor argument a window did not pass fell back to its constructor default, so opening a task and
    /// pressing Save changed the search without a word:</para>
    /// <list type="bullet">
    /// <item><description>every window reset <see cref="DigestionParams.RespectCleavageBlockingModifications"/> to false
    /// (no window has a control for it);</description></item>
    /// <item><description>the crosslink window reset <see cref="DigestionParams.MaxModsForPeptide"/> to 2 (it has no control
    /// for it);</description></item>
    /// <item><description>the calibration window showed the loaded <see cref="DigestionParams.InitiatorMethionineBehavior"/>
    /// in its combo box but saved Variable whatever the box said;</description></item>
    /// <item><description>the calibration, GPTMD, glyco and crosslink windows reset
    /// <see cref="DigestionParams.GeneratehUnlabeledProteinsForSilac"/> to true;</description></item>
    /// <item><description>every window reset KeepNGlycopeptide and KeepOGlycopeptide to false (nothing reads them today,
    /// but they are written to the settings file).</description></item>
    /// </list>
    /// <para><b>What the windows do now.</b> Each passes <see cref="TaskWindowDigestionParams.ForSave"/> only the values its
    /// controls show, and every other setting comes from the loaded task. The WPF windows cannot be referenced from this test
    /// project, so <see cref="WindowSaves"/> repeats each window's call. <b>Keep the two in step</b>: a window that starts
    /// showing a setting passes it in both places.</para>
    /// <para><see cref="ASettingNoWindowShows_SurvivesASave"/> walks the <see cref="DigestionParams"/> constructor by
    /// reflection, so a setting mzLib adds later fails there until <see cref="TaskWindowDigestionParams.ForSave"/> carries
    /// it. A hand-written list would not catch it: a forgotten setting falls back to its default, and a test that leaves the
    /// setting at its default agrees by accident. That is how RespectCleavageBlockingModifications got past
    /// FileSpecificDigestionParamsTests, which guards the same kind of rebuild in MetaMorpheusTask.</para>
    /// </remarks>
    [TestFixture]
    public static class TaskWindowDigestionParamsTests
    {
        /// <summary>The five settings every task window shows, and so always passes to <see cref="TaskWindowDigestionParams.ForSave"/>.</summary>
        private static readonly string[] ShownByEveryWindow =
            { "protease", "maxMissedCleavages", "minPeptideLength", "maxPeptideLength", "maxModificationIsoforms" };

        private static ConstructorInfo FullConstructor =>
            typeof(DigestionParams).GetConstructors().OrderByDescending(c => c.GetParameters().Length).First();

        /// <summary>Every <see cref="DigestionParams"/> constructor parameter that some window does not show.</summary>
        private static IEnumerable<string> SettingsSomeWindowDoesNotShow() =>
            FullConstructor.GetParameters().Select(p => p.Name).Where(name => !ShownByEveryWindow.Contains(name));

        /// <summary>
        /// A task whose only non-default setting is the named one, opened in a window that shows none of the optional settings
        /// and saved unchanged, keeps that setting. This covers every constructor parameter, including any mzLib adds later.
        /// </summary>
        [Test]
        [TestCaseSource(nameof(SettingsSomeWindowDoesNotShow))]
        public static void ASettingNoWindowShows_SurvivesASave(string setting)
        {
            ParameterInfo[] parameters = FullConstructor.GetParameters();
            object[] arguments = parameters.Select(p => p.Name == setting ? NotTheDefault(p) : p.DefaultValue).ToArray();
            var loaded = (DigestionParams)FullConstructor.Invoke(arguments);

            DigestionParams saved = TaskWindowDigestionParams.ForSave(loaded, TaskWindowSearchMode.ProteaseToShow(loaded).Name,
                loaded.MaxMissedCleavages, loaded.MinLength, loaded.MaxLength, loaded.MaxModificationIsoforms);

            AssertSameSettings(saved, loaded, $"{setting} was reset by saving a task in a window that does not show it");
        }

        /// <summary>
        /// Each window's call to <see cref="TaskWindowDigestionParams.ForSave"/>, with the window's controls showing what the
        /// loaded task holds (the user opened the task and changed nothing). Mirrors SaveButton_Click in each window.
        /// </summary>
        private static readonly Dictionary<string, Func<DigestionParams, DigestionParams>> WindowSaves = new()
        {
            ["Calibrate"] = loaded => TaskWindowDigestionParams.ForSave(loaded, loaded.SpecificProtease.Name,
                loaded.MaxMissedCleavages, loaded.MinLength, loaded.MaxLength, loaded.MaxModificationIsoforms,
                maxModsForPeptides: loaded.MaxMods,
                initiatorMethionineBehavior: loaded.InitiatorMethionineBehavior),

            ["GPTMD"] = loaded => TaskWindowDigestionParams.ForSave(loaded, loaded.SpecificProtease.Name,
                loaded.MaxMissedCleavages, loaded.MinLength, loaded.MaxLength, loaded.MaxModificationIsoforms,
                maxModsForPeptides: loaded.MaxMods,
                initiatorMethionineBehavior: loaded.InitiatorMethionineBehavior),

            ["Glyco"] = loaded => TaskWindowDigestionParams.ForSave(loaded, loaded.SpecificProtease.Name,
                loaded.MaxMissedCleavages, loaded.MinLength, loaded.MaxLength, loaded.MaxModificationIsoforms,
                maxModsForPeptides: loaded.MaxMods,
                initiatorMethionineBehavior: loaded.InitiatorMethionineBehavior,
                searchMode: TaskWindowSearchMode.ForSemiSpecificChoice(TaskWindowSearchMode.IsSemiSpecific(loaded))),

            // the crosslink window has no control for MaxModsForPeptide
            ["Crosslink"] = loaded => TaskWindowDigestionParams.ForSave(loaded, loaded.SpecificProtease.Name,
                loaded.MaxMissedCleavages, loaded.MinLength, loaded.MaxLength, loaded.MaxModificationIsoforms,
                initiatorMethionineBehavior: loaded.InitiatorMethionineBehavior),

            ["Search"] = loaded => TaskWindowDigestionParams.ForSave(loaded, loaded.SpecificProtease.Name,
                loaded.MaxMissedCleavages, loaded.MinLength, loaded.MaxLength, loaded.MaxModificationIsoforms,
                maxModsForPeptides: loaded.MaxMods,
                initiatorMethionineBehavior: loaded.InitiatorMethionineBehavior,
                searchMode: (loaded.SearchModeType, loaded.FragmentationTerminus),
                generateUnlabeledProteinsForSilac: loaded.GeneratehUnlabeledProteinsForSilac),
        };

        private static IEnumerable<string> Windows() => WindowSaves.Keys;

        /// <summary>
        /// Opening a task in any window and saving it without touching anything gives back the same digestion settings. The
        /// loaded task has every setting away from its default, so a setting the save forgets cannot agree by accident.
        /// </summary>
        [Test]
        [TestCaseSource(nameof(Windows))]
        public static void OpeningAndSavingATask_KeepsEveryDigestionSetting(string window)
        {
            Func<DigestionParams, DigestionParams> save = WindowSaves[window];
            var loaded = new DigestionParams("trypsin", maxMissedCleavages: 3, minPeptideLength: 8, maxPeptideLength: 40,
                maxModificationIsoforms: 512, initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave,
                maxModsForPeptides: 5, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.Both,
                generateUnlabeledProteinsForSilac: false, keepNGlycopeptide: true, keepOGlycopeptide: true,
                respectCleavageBlockingModifications: true);

            AssertSameSettings(save(loaded), loaded, $"the {window} window changed a setting the user did not touch");
        }

        /// <summary>
        /// What a window shows is what it saves: a value passed to <see cref="TaskWindowDigestionParams.ForSave"/> replaces the
        /// loaded task's value, so keeping hidden settings never undoes a change the user made.
        /// </summary>
        [Test]
        public static void AValueTheWindowShows_ReplacesTheLoadedOne()
        {
            var loaded = new DigestionParams("trypsin", initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave,
                maxModsForPeptides: 5, searchModeType: CleavageSpecificity.Semi, generateUnlabeledProteinsForSilac: false);

            DigestionParams saved = TaskWindowDigestionParams.ForSave(loaded, "Arg-C", 1, 6, 30, 256,
                maxModsForPeptides: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                searchMode: (CleavageSpecificity.Full, FragmentationTerminus.Both),
                generateUnlabeledProteinsForSilac: true);

            Assert.That(saved.SpecificProtease.Name, Is.EqualTo("Arg-C"));
            Assert.That(saved.MaxMissedCleavages, Is.EqualTo(1));
            Assert.That(saved.MinLength, Is.EqualTo(6));
            Assert.That(saved.MaxLength, Is.EqualTo(30));
            Assert.That(saved.MaxModificationIsoforms, Is.EqualTo(256));
            Assert.That(saved.MaxMods, Is.EqualTo(1));
            Assert.That(saved.InitiatorMethionineBehavior, Is.EqualTo(InitiatorMethionineBehavior.Retain));
            Assert.That(saved.SearchModeType, Is.EqualTo(CleavageSpecificity.Full));
            Assert.That(saved.GeneratehUnlabeledProteinsForSilac, Is.True);
        }

        /// <summary>
        /// A new task, or one without protein digestion parameters (an RNA task), has nothing to keep: the window saves the
        /// constructor defaults for every setting it does not show, exactly as before this fix.
        /// </summary>
        [Test]
        public static void ATaskWithNothingToKeep_SavesTheDefaults()
        {
            foreach (IDigestionParams loaded in new IDigestionParams[] { null, new RnaDigestionParams() })
            {
                DigestionParams saved = TaskWindowDigestionParams.ForSave(loaded, "trypsin", 2, 7, int.MaxValue, 1024);
                AssertSameSettings(saved, new DigestionParams("trypsin"), $"loaded {loaded?.GetType().Name ?? "nothing"}");
            }
        }

        /// <summary>A value of the parameter's type that is not its default, to show whether a save carries it.</summary>
        private static object NotTheDefault(ParameterInfo parameter)
        {
            Type type = parameter.ParameterType;
            object defaultValue = parameter.DefaultValue;
            if (type == typeof(bool))
                return !(bool)defaultValue;
            if (type == typeof(int))
                return (int)defaultValue == int.MaxValue ? 40 : (int)defaultValue + 1;
            if (type.IsEnum)
                return Enum.GetValues(type).Cast<object>().First(value => Convert.ToInt64(value) != Convert.ToInt64(defaultValue));

            Assert.Fail($"DigestionParams has a new constructor parameter, {parameter.Name}, of type {type.Name}. Teach " +
                        $"{nameof(NotTheDefault)} a non-default value for it, and make TaskWindowDigestionParams.ForSave carry it.");
            return null;
        }

        /// <summary>
        /// Compares every public property, so a failure names the setting that changed, then the whole object (Equals is what
        /// the index cache and the rest of MetaMorpheus use).
        /// </summary>
        private static void AssertSameSettings(DigestionParams saved, DigestionParams expected, string because)
        {
            foreach (PropertyInfo property in typeof(DigestionParams).GetProperties(BindingFlags.Public | BindingFlags.Instance))
            {
                Assert.That(property.GetValue(saved), Is.EqualTo(property.GetValue(expected)), $"{property.Name}: {because}");
            }
            Assert.That(saved, Is.EqualTo(expected), because);
        }
    }
}
