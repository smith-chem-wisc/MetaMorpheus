using System;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Reflection;
using GuiFunctions;
using NUnit.Framework;
using EngineLayer;
using MzLibUtil;
using GuiFunctions.Util;

namespace Test.GuiTests
{
    [TestFixture]
    public class GuiGlobalParamsTests
    {
        private static string SettingsPath => (string)typeof(GuiGlobalParamsViewModel)
            .GetProperty("SettingsPath", BindingFlags.NonPublic | BindingFlags.Static)
            .GetValue(null);

        [SetUp]
        public void SetUp()
        {
            // Ensure singleton reset before each test
            ResetSingleton();
            // Ensure GlobalVariables initialized for tests
            GlobalVariables.SetUpGlobalVariables();
            
            // Clear any event handlers from previous tests
            ClearEventHandlers();
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
            ClearEventHandlers();
        }

        private static void ResetSingleton()
        {
            var vmType = typeof(GuiGlobalParamsViewModel);
            var instField = vmType.GetField("_instance", BindingFlags.NonPublic | BindingFlags.Static);
            instField.SetValue(null, null);
        }

        private static void ClearEventHandlers()
        {
            // Clear static event handlers to prevent cross-test contamination
            var requestModeSwitchField = typeof(GuiGlobalParamsViewModel)
                .GetProperty("RequestModeSwitchConfirmation", BindingFlags.Public | BindingFlags.Static);
            requestModeSwitchField.SetValue(null, null);
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
                // Skip IsRnaMode as it has special behavior tested separately
                if (prop.Name == nameof(GuiGlobalParamsViewModel.IsRnaMode)) continue;

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
            GlobalVariables.SetUpGlobalVariables();

            // Hook up event handler that auto-approves the switch
            EventHandler<ModeSwitchRequestEventArgs> handler = (s, e) =>
            {
                e.Result = ModeSwitchResult.SwitchKeepFiles;
            };
            GuiGlobalParamsViewModel.RequestModeSwitchConfirmation += handler;

            int eventCount = 0;
            vm.PropertyChanged += (s, e) => 
            { 
                if (e.PropertyName == nameof(vm.IsRnaMode) || e.PropertyName == nameof(vm.MainWindowTitle)) 
                    eventCount++; 
            };

            try
            {
                vm.IsRnaMode = true;
                Assert.That(GlobalVariables.AnalyteType, Is.EqualTo(EngineLayer.AnalyteType.Oligo));
                Assert.That(vm.MainWindowTitle, Does.Contain("RNA"));

                vm.IsRnaMode = false;
                Assert.That(GlobalVariables.AnalyteType, Is.EqualTo(EngineLayer.AnalyteType.Peptide));
                Assert.That(vm.MainWindowTitle, Does.Contain("Protein"));

                Assert.That(eventCount, Is.GreaterThanOrEqualTo(2));
            }
            finally
            {
                GuiGlobalParamsViewModel.RequestModeSwitchConfirmation -= handler;
            }
        }

        [Test]
        public void IsRnaMode_WithCancel_DoesNotSwitch()
        {
            var vm = GuiGlobalParamsViewModel.Instance;
            GlobalVariables.SetUpGlobalVariables();

            // Hook up event handler that cancels the switch
            EventHandler<ModeSwitchRequestEventArgs> handler = (s, e) =>
            {
                e.Result = ModeSwitchResult.Cancel;
            };
            GuiGlobalParamsViewModel.RequestModeSwitchConfirmation += handler;

            try
            {
                bool originalMode = vm.IsRnaMode;
                var originalAnalyteType = GlobalVariables.AnalyteType;

                // Try to switch
                vm.IsRnaMode = !originalMode;

                // Should not have changed
                Assert.That(vm.IsRnaMode, Is.EqualTo(originalMode));
                Assert.That(GlobalVariables.AnalyteType, Is.EqualTo(originalAnalyteType));
            }
            finally
            {
                GuiGlobalParamsViewModel.RequestModeSwitchConfirmation -= handler;
            }
        }

        [Test]
        public void IsRnaMode_WithRememberDecision_UpdatesSettings()
        {
            var vm = GuiGlobalParamsViewModel.Instance;
            GlobalVariables.SetUpGlobalVariables();

            Assert.That(vm.AskAboutModeSwitch, Is.True, "Should start with AskAboutModeSwitch = true");

            // Hook up event handler that approves and remembers
            EventHandler<ModeSwitchRequestEventArgs> handler = (s, e) =>
            {
                e.Result = ModeSwitchResult.SwitchRemoveFiles;
                e.RememberMyDecision = true;
            };
            GuiGlobalParamsViewModel.RequestModeSwitchConfirmation += handler;

            try
            {
                vm.IsRnaMode = true;

                // Should have updated the cached result and disabled asking
                Assert.That(vm.AskAboutModeSwitch, Is.False);
                Assert.That(vm.CachedModeSwitchResult, Is.EqualTo(ModeSwitchResult.SwitchRemoveFiles));
            }
            finally
            {
                GuiGlobalParamsViewModel.RequestModeSwitchConfirmation -= handler;
            }
        }

        [Test]
        public void IsRnaMode_EventInvoked_WithCorrectArguments()
        {
            var vm = GuiGlobalParamsViewModel.Instance;
            GlobalVariables.SetUpGlobalVariables();

            bool eventWasInvoked = false;
            ModeSwitchRequestEventArgs capturedArgs = null;

            EventHandler<ModeSwitchRequestEventArgs> handler = (s, e) =>
            {
                eventWasInvoked = true;
                capturedArgs = e;
                e.Result = ModeSwitchResult.SwitchKeepFiles;
            };
            GuiGlobalParamsViewModel.RequestModeSwitchConfirmation += handler;

            try
            {
                vm.IsRnaMode = true;

                Assert.That(eventWasInvoked, Is.True, "Event should have been invoked");
                Assert.That(capturedArgs, Is.Not.Null, "Event args should not be null");
                Assert.That(capturedArgs.Result, Is.EqualTo(ModeSwitchResult.SwitchKeepFiles));
            }
            finally
            {
                GuiGlobalParamsViewModel.RequestModeSwitchConfirmation -= handler;
            }
        }

        [Test]
        public void CachedModeSwitchResult_GetSet_WorksCorrectly()
        {
            var vm = GuiGlobalParamsViewModel.Instance;

            // Test all enum values
            foreach (ModeSwitchResult result in Enum.GetValues(typeof(ModeSwitchResult)))
            {
                vm.CachedModeSwitchResult = result;
                Assert.That(vm.CachedModeSwitchResult, Is.EqualTo(result));
            }
        }

        [Test]
        public void CachedModeSwitchResult_TriggersPropertyChanged()
        {
            var vm = GuiGlobalParamsViewModel.Instance;

            bool eventFired = false;
            PropertyChangedEventHandler handler = (s, e) =>
            {
                if (e.PropertyName == nameof(vm.CachedModeSwitchResult))
                    eventFired = true;
            };

            vm.PropertyChanged += handler;
            try
            {
                vm.CachedModeSwitchResult = ModeSwitchResult.SwitchKeepFiles;
                Assert.That(eventFired, Is.True);
            }
            finally
            {
                vm.PropertyChanged -= handler;
            }
        }

        [Test]
        public void AskAboutModeSwitch_GetSet_WorksCorrectly()
        {
            var vm = GuiGlobalParamsViewModel.Instance;

            vm.AskAboutModeSwitch = false;
            Assert.That(vm.AskAboutModeSwitch, Is.False);

            vm.AskAboutModeSwitch = true;
            Assert.That(vm.AskAboutModeSwitch, Is.True);
        }

        [Test]
        public void AskAboutModeSwitch_TriggersPropertyChanged()
        {
            var vm = GuiGlobalParamsViewModel.Instance;

            bool eventFired = false;
            PropertyChangedEventHandler handler = (s, e) =>
            {
                if (e.PropertyName == nameof(vm.AskAboutModeSwitch))
                    eventFired = true;
            };

            vm.PropertyChanged += handler;
            try
            {
                vm.AskAboutModeSwitch = !vm.AskAboutModeSwitch;
                Assert.That(eventFired, Is.True);
            }
            finally
            {
                vm.PropertyChanged -= handler;
            }
        }

        [Test]
        public void AllModeSwitchValues_ContainsAllEnumValues()
        {
            var vm = GuiGlobalParamsViewModel.Instance;
            var allValues = Enum.GetValues(typeof(ModeSwitchResult)).Cast<ModeSwitchResult>().ToList();

            Assert.That(vm.AllModeSwitchValues.Count, Is.EqualTo(allValues.Count));
            
            foreach (var value in allValues)
            {
                Assert.That(vm.AllModeSwitchValues.Contains(value), Is.True, 
                    $"AllModeSwitchValues should contain {value}");
            }
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

        [Test]
        public void GuiGlobalParams_Equals_ChecksAllModeSwitchFields()
        {
            var a = new GuiGlobalParams();
            var b = a.Clone();

            Assert.That(a.Equals(b), Is.True, "Cloned params should be equal");

            // Change IsRnaMode
            b.IsRnaMode = !a.IsRnaMode;
            Assert.That(a.Equals(b), Is.False, "Should not be equal after changing IsRnaMode");
            b.IsRnaMode = a.IsRnaMode;

            // Change AskAboutModeSwitch
            b.AskAboutModeSwitch = !a.AskAboutModeSwitch;
            Assert.That(a.Equals(b), Is.False, "Should not be equal after changing AskAboutModeSwitch");
            b.AskAboutModeSwitch = a.AskAboutModeSwitch;

            // Change CachedModeSwitchResult
            b.CachedModeSwitchResult = ModeSwitchResult.SwitchRemoveFiles;
            Assert.That(a.Equals(b), Is.False, "Should not be equal after changing CachedModeSwitchResult");
            b.CachedModeSwitchResult = a.CachedModeSwitchResult;

            Assert.That(a.Equals(b), Is.True, "Should be equal after reverting all changes");
        }

        [Test]
        public void GuiGlobalParams_Clone_CopiesAllModeSwitchFields()
        {
            var original = new GuiGlobalParams
            {
                IsRnaMode = true,
                AskAboutModeSwitch = false,
                CachedModeSwitchResult = ModeSwitchResult.SwitchRemoveFiles
            };

            var clone = original.Clone();

            Assert.That(clone.IsRnaMode, Is.EqualTo(original.IsRnaMode));
            Assert.That(clone.AskAboutModeSwitch, Is.EqualTo(original.AskAboutModeSwitch));
            Assert.That(clone.CachedModeSwitchResult, Is.EqualTo(original.CachedModeSwitchResult));
        }

        [Test]
        public void ModeSwitchRequestEventArgs_DefaultValues()
        {
            var args = new ModeSwitchRequestEventArgs();

            Assert.That(args.Result, Is.EqualTo(ModeSwitchResult.Cancel), 
                "Default result should be Cancel");
            Assert.That(args.RememberMyDecision, Is.False, 
                "Default RememberMyDecision should be false");
        }

        [Test]
        public void ModeSwitchRequestEventArgs_CanSetValues()
        {
            var args = new ModeSwitchRequestEventArgs
            {
                Result = ModeSwitchResult.SwitchRemoveFiles,
                RememberMyDecision = true
            };

            Assert.That(args.Result, Is.EqualTo(ModeSwitchResult.SwitchRemoveFiles));
            Assert.That(args.RememberMyDecision, Is.True);
        }

        [Test]
        public void IsRnaMode_SavesAndLoadsCorrectly()
        {
            var vm = GuiGlobalParamsViewModel.Instance;

            // Setup event handler
            EventHandler<ModeSwitchRequestEventArgs> handler = (s, e) =>
            {
                e.Result = ModeSwitchResult.SwitchKeepFiles;
            };
            GuiGlobalParamsViewModel.RequestModeSwitchConfirmation += handler;

            try
            {
                // Set RNA mode
                vm.IsRnaMode = true;
                vm.Save();

                // Reset singleton and reload
                ResetSingleton();
                var vm2 = GuiGlobalParamsViewModel.Instance;

                Assert.That(vm2.IsRnaMode, Is.True, "IsRnaMode should persist after save/load");
            }
            finally
            {
                GuiGlobalParamsViewModel.RequestModeSwitchConfirmation -= handler;
            }
        }

        [Test]
        public void ModeSwitchSettings_SaveAndLoad_PersistCorrectly()
        {
            var vm = GuiGlobalParamsViewModel.Instance;

            // Set mode switch settings
            vm.AskAboutModeSwitch = false;
            vm.CachedModeSwitchResult = ModeSwitchResult.SwitchRemoveFiles;
            vm.Save();

            // Reset and reload
            ResetSingleton();
            var vm2 = GuiGlobalParamsViewModel.Instance;

            Assert.That(vm2.AskAboutModeSwitch, Is.False);
            Assert.That(vm2.CachedModeSwitchResult, Is.EqualTo(ModeSwitchResult.SwitchRemoveFiles));
        }

        [Test]
        public void ModeSwitchRequestEventArgs_IsInheritedFromEventArgs()
        {
            var args = new ModeSwitchRequestEventArgs();
            Assert.That(args, Is.InstanceOf<EventArgs>());
        }

        [Test]
        public void ModeSwitchResult_EnumValues_AreCorrect()
        {
            var values = Enum.GetValues(typeof(ModeSwitchResult)).Cast<ModeSwitchResult>().ToList();
            
            Assert.That(values, Contains.Item(ModeSwitchResult.Cancel));
            Assert.That(values, Contains.Item(ModeSwitchResult.SwitchKeepFiles));
            Assert.That(values, Contains.Item(ModeSwitchResult.SwitchRemoveFiles));
            Assert.That(values.Count, Is.EqualTo(3), "ModeSwitchResult should have exactly 3 values");
        }
    }
}
