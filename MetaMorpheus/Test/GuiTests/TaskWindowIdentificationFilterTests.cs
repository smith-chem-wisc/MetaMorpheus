using EngineLayer;
using EngineLayer.SpectrumMatch;
using GuiFunctions.Util;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading;
using TaskLayer;

namespace Test.GuiTests
{
    /// <summary>
    /// What the task windows show and save for the identification filter and its threshold: the filter choice and one
    /// threshold box.
    /// </summary>
    /// <remarks>
    /// The Search window loads them with <see cref="TaskWindowIdentificationFilter.Load(CommonParameters, bool)"/> and saves
    /// them with <see cref="TaskWindowIdentificationFilter.ForSave"/>. The WPF windows cannot be referenced from this test
    /// project, so <see cref="WindowSaves"/> repeats each window's calls. <b>Keep the two in step</b>: a window that starts
    /// showing the filter is added here.
    /// </remarks>
    [TestFixture]
    public static class TaskWindowIdentificationFilterTests
    {
        /// <summary>
        /// Each window's load and save of the filter, with its controls showing what the loaded task holds (the user opened the
        /// task and changed nothing), building the CommonParameters it saves. Mirrors UpdateFieldsFromTask and SaveButton_Click;
        /// the window passes the RNA mode the GUI is in.
        /// </summary>
        private static readonly Dictionary<string, Func<CommonParameters, bool, CommonParameters>> WindowSaves = new()
        {
            ["Search"] = OpenAndSave,
        };

        private static CommonParameters OpenAndSave(CommonParameters loaded, bool isRnaMode)
        {
            var (shownFilter, shownThreshold) = TaskWindowIdentificationFilter.Load(loaded, isRnaMode);
            var (filter, threshold) = TaskWindowIdentificationFilter.ForSave(shownFilter, shownThreshold);
            return new CommonParameters(digestionParams: loaded.DigestionParams, identificationFilter: filter, filterThreshold: threshold);
        }

        private static IEnumerable<string> Windows() => WindowSaves.Keys;

        private static IEnumerable<FilterType> Filters() => Enum.GetValues<FilterType>();

        private static IEnumerable<TestCaseData> WindowsAndFilters() =>
            Windows().SelectMany(window => Filters().Select(filter => new TestCaseData(window, filter)));

        /// <summary>Every filter the enum holds is offered by the windows, once, with a name and what it means.</summary>
        [Test]
        [TestCaseSource(nameof(Filters))]
        public static void EveryFilter_IsOfferedOnce_WithANameAndWhatItMeans(FilterType filter)
        {
            var offered = TaskWindowIdentificationFilter.Choices.Where(c => c.Filter == filter).ToList();
            Assert.That(offered, Has.Count.EqualTo(1));
            Assert.That(offered[0].Label, Is.Not.Empty);
            Assert.That(offered[0].Description, Is.Not.Empty);
        }

        [Test]
        public static void TheChoices_AreNamedAsTheMaintainerAskedForThem()
        {
            Assert.That(TaskWindowIdentificationFilter.Choices.Select(c => (c.Filter, c.Label)), Is.EqualTo(new[]
            {
                (FilterType.PepQValue, "PEP q-value"),
                (FilterType.QValueNotch, "q-value notch"),
                (FilterType.QValue, "q-value (pooled)"),
            }));
            Assert.That(TaskWindowIdentificationFilter.Choices.Single(c => c.Filter == FilterType.PepQValue).Description,
                Does.Contain("q-value notch"), "the PEP q-value says it falls back to the notch where PEP could not be trained");
        }

        /// <summary>
        /// Opening a task and saving it without touching anything keeps its filter and threshold, for every filter. The
        /// threshold is away from its default, so a window that saved the default could not agree by accident.
        /// </summary>
        [Test]
        [TestCaseSource(nameof(WindowsAndFilters))]
        public static void OpeningAndSavingATask_KeepsItsFilterAndThreshold(string window, FilterType filter)
        {
            var loaded = new CommonParameters(identificationFilter: filter, filterThreshold: 0.02);

            CommonParameters saved = WindowSaves[window](loaded, false);

            Assert.That((saved.IdentificationFilter, saved.FilterThreshold), Is.EqualTo((filter, 0.02)),
                $"the {window} window changed a filter the user did not touch");
        }

        /// <summary>The q-value notch task at the old default, 0.01, is not moved to the new bottom-up default by a save.</summary>
        [Test]
        [TestCaseSource(nameof(Windows))]
        public static void ANotchTaskAt001_StaysTheNotch(string window)
        {
            var loaded = new CommonParameters(identificationFilter: FilterType.QValueNotch);

            CommonParameters saved = WindowSaves[window](loaded, false);

            Assert.That((saved.IdentificationFilter, saved.FilterThreshold), Is.EqualTo((FilterType.QValueNotch, 0.01)));
        }

        /// <summary>An old GUI PEP task (QValueThreshold 1.0, PepQValueThreshold 0.01, no filter named) stays a PEP task.</summary>
        [Test]
        [TestCaseSource(nameof(Windows))]
        public static void AnOldGuiPepTask_OpenedAndSaved_StaysThePepFilter(string window)
        {
            var old = Nett.Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Task2-SearchTaskconfig.toml"),
                MetaMorpheusTask.tomlConfig).CommonParameters;

            CommonParameters saved = WindowSaves[window](old, false);

            Assert.That((saved.IdentificationFilter, saved.FilterThreshold), Is.EqualTo((FilterType.PepQValue, 0.01)));
        }

        /// <summary>
        /// A task still at the bottom-up default (the PEP q-value at 0.01) shows the default for its search type: the q-value
        /// notch with a top-down protease or in RNA mode. Any other setting is shown as it is.
        /// </summary>
        [Test]
        public static void ATaskAtTheBottomUpDefault_ShowsTheDefaultForItsSearchType()
        {
            var topDown = new DigestionParams(protease: "top-down");
            var atDefault = new CommonParameters(digestionParams: topDown, identificationFilter: FilterType.PepQValue, filterThreshold: 0.01);
            Assert.That(TaskWindowIdentificationFilter.Load(atDefault, isRna: false), Is.EqualTo((FilterType.QValueNotch, "0.01")));

            var bottomUpAtDefault = new CommonParameters(identificationFilter: FilterType.PepQValue, filterThreshold: 0.01);
            Assert.That(TaskWindowIdentificationFilter.Load(bottomUpAtDefault, isRna: false), Is.EqualTo((FilterType.PepQValue, "0.01")));
            Assert.That(TaskWindowIdentificationFilter.Load(bottomUpAtDefault, isRna: true), Is.EqualTo((FilterType.QValueNotch, "0.01")));

            var chosenThreshold = new CommonParameters(digestionParams: topDown, identificationFilter: FilterType.PepQValue, filterThreshold: 0.02);
            Assert.That(TaskWindowIdentificationFilter.Load(chosenThreshold, isRna: false), Is.EqualTo((FilterType.PepQValue, "0.02")));

            var oldDefaultNotch = new CommonParameters(identificationFilter: FilterType.QValueNotch);
            Assert.That(TaskWindowIdentificationFilter.Load(oldDefaultNotch, isRna: false), Is.EqualTo((FilterType.QValueNotch, "0.01")),
                "the notch at 0.01 is never moved to the PEP q-value");
        }

        /// <summary>The top-down recommended settings set the q-value notch at 0.01, and the prompt says so.</summary>
        [Test]
        [TestCase(false)]
        [TestCase(true)]
        public static void TheTopDownRecommendation_IsTheNotchAt001(bool isRna)
        {
            Assert.That(TaskWindowIdentificationFilter.TopDownRecommended(isRna), Is.EqualTo((FilterType.QValueNotch, "0.01")));
            Assert.That(TaskWindowIdentificationFilter.TopDownRecommended(isRna).Filter, Is.EqualTo(CommonParameters.DefaultIdentificationFilter(isRna, true)));
            Assert.That(TaskWindowIdentificationFilter.TopDownRecommendedSentence(isRna),
                Is.EqualTo("Use 'q-value notch' at '0.01' for 'Identification Filter'"));
        }

        /// <summary>The threshold is shown and read in the invariant culture, whatever the machine's culture.</summary>
        [Test]
        public static void TheThreshold_IsShownAndReadInTheInvariantCulture()
        {
            CultureInfo original = Thread.CurrentThread.CurrentCulture;
            try
            {
                Thread.CurrentThread.CurrentCulture = new CultureInfo("de-DE");
                var loaded = new CommonParameters(identificationFilter: FilterType.QValueNotch, filterThreshold: 0.005);
                Assert.That(TaskWindowIdentificationFilter.Load(loaded, isRna: false), Is.EqualTo((FilterType.QValueNotch, "0.005")));
                Assert.That(TaskWindowIdentificationFilter.ForSave(FilterType.QValueNotch, "0.005"), Is.EqualTo((FilterType.QValueNotch, 0.005)));
            }
            finally
            {
                Thread.CurrentThread.CurrentCulture = original;
            }
        }

        /// <summary>The one threshold is above 0 and at most 1; at 1, which keeps every match, the window warns.</summary>
        [Test]
        [TestCase("0.01", false, false)]
        [TestCase("0.005", false, false)]
        [TestCase("1e-3", false, false)]
        [TestCase("1", false, true)]
        [TestCase("1.0", false, true)]
        [TestCase("0", true, false)]
        [TestCase("-0.01", true, false)]
        [TestCase("1.5", true, false)]
        [TestCase("", true, false)]
        [TestCase("abc", true, false)]
        [TestCase("0,01", true, false)]
        public static void TheThreshold_IsChecked(string text, bool refused, bool warned)
        {
            Assert.That(TaskWindowIdentificationFilter.ThresholdError(text) != null, Is.EqualTo(refused));
            if (!refused)
            {
                Assert.That(TaskWindowIdentificationFilter.ThresholdWarning(text) != null, Is.EqualTo(warned));
            }
        }
    }
}
