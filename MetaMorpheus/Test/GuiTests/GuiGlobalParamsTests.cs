using System;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Reflection;
using GuiFunctions;
using NUnit.Framework;
using EngineLayer;
using MzLibUtil;

namespace Test.GuiTests
{
    [TestFixture]
    public class GuiGlobalParamsTests
    {
        private static string SettingsPath => (string)typeof(GuiGlobalParamsViewModel)
            .GetField("SettingsPath", BindingFlags.NonPublic | BindingFlags.Static)
            .GetValue(null);

        [SetUp]
        public void SetUp()
        {
            // Ensure singleton reset before each test
            ResetSingleton();
            // Ensure GlobalVariables initialized for tests
            GlobalVariables.SetUpGlobalVariables();
        }

        [TearDown]
        public void TearDown()
        {
            // cleanup created settings file to avoid cross-test interference
            try
            {
                if (File.Exists(SettingsPath))
                {
                    File.Delete(SettingsPath);
                }
            }
            catch { }

            ResetSingleton();
        }

        private static void ResetSingleton()
        {
            var vmType = typeof(GuiGlobalParamsViewModel);
            var instField = vmType.GetField("_instance", BindingFlags.NonPublic | BindingFlags.Static);
            instField.SetValue(null, null);
        }

        [Test]
        public void Load_CreatesDefaultSettingsFile_WhenMissing()
        {
            if (File.Exists(SettingsPath)) File.Delete(SettingsPath);
            Assert.That(File.Exists(SettingsPath), Is.False);

            var vm = GuiGlobalParamsViewModel.Instance; // triggers Load()

            Assert.That(File.Exists(SettingsPath), Is.True, "Load should create settings file when missing");

            // ProteomeDirectory should be set to default when file created
            var expectedDefault = Path.Combine(GlobalVariables.DataDir, "Proteomes");
            Assert.That(vm.ProteomeDirectory, Is.EqualTo(expectedDefault));
        }

        [Test]
        public void Load_ReadsExistingTomlFile()
        {
            // create a TOML settings file with custom values
            string toml = "UserSpecifiedProteomeDir = \"/tmp/proteomes_custom\"\nAskAboutUpdating = false\nDecoyIdentifier = \"XDECOY\"\n";
            File.WriteAllText(SettingsPath, toml);

            var vm = GuiGlobalParamsViewModel.Instance; // Load reads file

            Assert.That(vm.ProteomeDirectory, Is.EqualTo("/tmp/proteomes_custom"));
            Assert.That(vm.AskAboutUpdating, Is.False);
            Assert.That(GlobalVariables.DecoyIdentifier, Is.EqualTo("XDECOY"));
        }

        [Test]
        public void Load_WithMalformedToml_FallsBackToDefaults()
        {
            // write malformed content that should cause Toml.ReadFile to throw
            File.WriteAllText(SettingsPath, "this is not toml = [}{");

            // Reset and instantiate
            ResetSingleton();
            var vm = GuiGlobalParamsViewModel.Instance; // should catch and rewrite defaults

            // Should have default proteome dir
            var expectedDefault = Path.Combine(GlobalVariables.DataDir, "Proteomes");
            Assert.That(vm.ProteomeDirectory, Is.EqualTo(expectedDefault));

            // File should exist (rewritten)
            Assert.That(File.Exists(SettingsPath), Is.True);
        }

        [Test]
        public void Save_WritesFile_And_IsDirtyBehavior()
        {
            // ensure instance created and loaded
            var vm = GuiGlobalParamsViewModel.Instance;
            Assert.That(vm.IsDirty(), Is.False);

            // toggle a property
            bool original = vm.AskAboutUpdating;
            vm.AskAboutUpdating = !original;
            Assert.That(vm.IsDirty(), Is.True, "After changing a property, IsDirty should be true");

            vm.Save();
            Assert.That(vm.IsDirty(), Is.False, "After Save, IsDirty should be false");

            // cleanup
            File.Delete(SettingsPath);
        }

        [Test]
        public void SettingsFileExists_ReturnsTrueAfterSave()
        {
            var vm = GuiGlobalParamsViewModel.Instance;
            vm.Save();
            Assert.That(GuiGlobalParamsViewModel.SettingsFileExists(), Is.True);
            File.Delete(SettingsPath);
        }

        [Test]
        public void ProteomeDirectory_GetterSetsDefaultWhenNull()
        {
            var vm = GuiGlobalParamsViewModel.Instance;
            // set _current.ProteomeDirectory to null via reflection to test getter branch
            var vmType = typeof(GuiGlobalParamsViewModel);
            var currentField = vmType.GetField("_current", BindingFlags.NonPublic | BindingFlags.Instance);
            var current = currentField.GetValue(vm);
            var backingField = current.GetType().GetField("<ProteomeDirectory>k__BackingField", BindingFlags.NonPublic | BindingFlags.Instance);
            backingField.SetValue(current, null);

            var expectedDefault = Path.Combine(GlobalVariables.DataDir, "Proteomes");
            Assert.That(vm.ProteomeDirectory, Is.EqualTo(expectedDefault));
        }

        [Test]
        public void ProteomeDirectory_Setter_OnlyUpdatesWhenDirectoryExists()
        {
            var vm = GuiGlobalParamsViewModel.Instance;
            string tmp = Path.Combine(TestContext.CurrentContext.TestDirectory, "tmpProteomes");
            if (Directory.Exists(tmp)) Directory.Delete(tmp, true);

            // setting to non-existing path should not change
            string before = vm.ProteomeDirectory;
            vm.ProteomeDirectory = tmp; // does not exist
            Assert.That(vm.ProteomeDirectory, Is.EqualTo(before));

            // setting to null path should not change
            vm.ProteomeDirectory = null;
            Assert.That(vm.ProteomeDirectory, Is.EqualTo(before));

            // create directory and set
            Directory.CreateDirectory(tmp);
            try
            {
                vm.ProteomeDirectory = tmp;
                Assert.That(vm.ProteomeDirectory, Is.EqualTo(tmp));
            }
            finally
            {
                Directory.Delete(tmp, true);
            }
        }

        [Test]
        public void BooleanProperties_SetAndGet_TriggerPropertyChanged()
        {
            var vm = GuiGlobalParamsViewModel.Instance;
            var vmType = typeof(GuiGlobalParamsViewModel);
            var props = vmType.GetProperties(BindingFlags.Public | BindingFlags.Instance);

            int testedCount = 0;

            foreach (var prop in props)
            {
                // skip non-bool properties, non-read/write, and indexers
                if (prop.PropertyType != typeof(bool)) continue;
                if (!prop.CanRead || !prop.CanWrite) continue;
                if (prop.GetIndexParameters().Length > 0) continue;

                testedCount++;

                bool eventFired = false;
                PropertyChangedEventHandler handler = (s, e) =>
                {
                    if (e.PropertyName == prop.Name) eventFired = true;
                };

                vm.PropertyChanged += handler;
                try
                {
                    // read original value
                    bool before = (bool)prop.GetValue(vm);
                    // set to opposite
                    prop.SetValue(vm, !before);
                    // assert value changed
                    var afterObj = prop.GetValue(vm);
                    Assert.That(afterObj, Is.EqualTo((object)!before), $"Property {prop.Name} did not toggle as expected.");
                    // assert event fired
                    Assert.That(eventFired, Is.True, $"Changing {prop.Name} should raise PropertyChanged for that property.");
                    // restore original to avoid side-effects
                    prop.SetValue(vm, before);
                }
                finally
                {
                    vm.PropertyChanged -= handler;
                }
            }

            Assert.That(testedCount, Is.GreaterThan(0), "No boolean properties were found to test on GuiGlobalParamsViewModel.");
        }

        [Test]
        public void IsRnaMode_TogglesGlobalAnalyteType_And_MainWindowTitle_And_RaisesEvents()
        {
            var vm = GuiGlobalParamsViewModel.Instance;
            // capture previous analyte
            GlobalVariables.SetUpGlobalVariables(); // ensure defaults

            int eventCount = 0;
            vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.IsRnaMode) || e.PropertyName == nameof(vm.MainWindowTitle)) eventCount++; };

            vm.IsRnaMode = true;
            Assert.That(GlobalVariables.AnalyteType, Is.EqualTo(EngineLayer.AnalyteType.Oligo));
            Assert.That(vm.MainWindowTitle, Does.Contain("RNA"));

            vm.IsRnaMode = false;
            Assert.That(GlobalVariables.AnalyteType, Is.EqualTo(EngineLayer.AnalyteType.Peptide));
            Assert.That(vm.MainWindowTitle, Does.Contain("Protein"));

            Assert.That(eventCount, Is.GreaterThanOrEqualTo(2));
        }

        [Test]
        public void GuiGlobalParams_Equals_And_Clone_WorkAsExpected()
        {
            var a = new GuiGlobalParams();
            var b = a.Clone();
            Assert.That(a.Equals(b), Is.True);

            // change a field on b by writing to its backing field
            var field = typeof(GuiGlobalParams).GetField("<UseTopDownParams>k__BackingField", BindingFlags.NonPublic | BindingFlags.Instance);
            bool original = (bool)field.GetValue(b);
            field.SetValue(b, !original);

            Assert.That(a.Equals(b), Is.False);
            Assert.That(a.Equals(null), Is.False);
        }
    }
}
