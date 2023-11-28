using EngineLayer;
using GuiFunctions;
using Nett;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using TaskLayer;

namespace Test
{

    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SaveSettingsTests
    {
        static SaveSettingsTests()
        {
            GlobalVariables.SetUpGlobalVariables();
        }
        public static IEnumerable<MetaMorpheusTask> SwitchTestCases =>
            new List<MetaMorpheusTask>
            {
                new SearchTask(),
                new SpectralAveragingTask(),
                new CalibrationTask(),
                new GlycoSearchTask(),
                new XLSearchTask(),
                new GptmdTask()
            };

        [OneTimeSetUp]
        public static void Setup() {}

        [TearDown]
        public static void TearDown()
        {
            String path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DefaultParameters");
            if (Directory.Exists(path))
            {
                Directory.Delete(path, true);
            }
        }

        /// <summary>
        /// fake window for the test
        /// </summary>
        private class FakeGuiWindow
        {
            public MetaMorpheusTask TheTask { get; internal set; }

            public FakeGuiWindow(MyTask taskType)
            {
                TheTask = taskType switch
                {
                    MyTask.Search => new SearchTask(),
                    MyTask.Gptmd => new GptmdTask(),
                    MyTask.Calibrate => new CalibrationTask(),
                    MyTask.XLSearch => new XLSearchTask(),
                    MyTask.GlycoSearch => new GlycoSearchTask(),
                    MyTask.Average => new SpectralAveragingTask(),
                    _ => throw new ArgumentOutOfRangeException(nameof(taskType), taskType, null)
                };
            }

            public MetaMorpheusTask GetTaskFromGui()
            {
                return TheTask;
            }

            public void UpdateFieldsFromNewTask(MetaMorpheusTask task)
            {
                TheTask = task;
            }
        }

        /// <summary>
        /// test the constructor without the defaultparameter directory
        /// </summary>
        /// <param name="task">task currenty working on</param>
        [Test]
        [NonParallelizable]
        [TestCaseSource(nameof(SwitchTestCases))]
        public void ConstructorTestWithoutDefaultParameters(MetaMorpheusTask task)
        {
            var fakeTaskWindow = new FakeGuiWindow(task.TaskType);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(task,
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var testDictionary = testTaskSettingViewModel.AllSettingsDict;
            String pathToCheck = Path.Combine(GlobalVariables.DataDir, "DefaultParameters");
            Assert.That(Directory.Exists(pathToCheck));
            Assert.That(testDictionary.ContainsKey("DefaultSetting(Default)"));
            var testWindow = fakeTaskWindow.GetTaskFromGui();
            var comparingWindow = testDictionary["DefaultSetting(Default)"];
            Assert.That(!testWindow.Equals(null));
            TomlFileFolderSerializer.Save("testFile", testWindow);
            string pathOfTestFile = TomlFileFolderSerializer.GetFilePath(testWindow.GetType(), "testFile");
            string pathOfComparingFile =
                TomlFileFolderSerializer.GetFilePath(comparingWindow.GetType(), "DefaultSetting(Default)");
            var linesA = File.ReadAllLines(pathOfTestFile);
            var linesB = File.ReadAllLines(pathOfComparingFile);

            Assert.That(linesA.Length, Is.EqualTo(linesB.Length));
            for (int i = 0; i < linesA.Length; i++)
            {
                Assert.That(linesA[i], Is.EqualTo(linesB[i]));
            }
        }

        /// <summary>
        /// test the constructor with the defaultparameter directory but nothing is in it
        /// </summary>
        /// <param name="task">task currenty working on</param>
        [Test]
        [NonParallelizable]
        public void ConstructorTestWithDefaultParametersWithoutDefaultSettings()
        {
            //setup
            Directory.CreateDirectory(Path.Combine(GlobalVariables.DataDir, "DefaultParameters"));
            //test
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(),
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var testDictionary = testTaskSettingViewModel.AllSettingsDict;
            Assert.That(testDictionary.ContainsKey("DefaultSetting(Default)"));
            var testWindow = fakeTaskWindow.GetTaskFromGui();
            var comparingWindow = testDictionary["DefaultSetting(Default)"];
            Assert.That(!testWindow.Equals(null));
            TomlFileFolderSerializer.Save("testFile", testWindow);
            string pathOfTestFile = TomlFileFolderSerializer.GetFilePath(testWindow.GetType(), "testFile");
            string pathOfComparingFile =
                TomlFileFolderSerializer.GetFilePath(comparingWindow.GetType(), "DefaultSetting(Default)");
            var linesA = File.ReadAllLines(pathOfTestFile);
            var linesB = File.ReadAllLines(pathOfComparingFile);

            Assert.That(linesA.Length, Is.EqualTo(linesB.Length));
            for (int i = 0; i < linesA.Length; i++)
            {
                Assert.That(linesA[i], Is.EqualTo(linesB[i]));
            }
        }

        /// <summary>
        /// test the constructor with the defaultparameter tasks
        /// </summary>
        /// <param name="task">task currenty working on</param>
        [Test]
        [NonParallelizable]
        [TestCaseSource(nameof(SwitchTestCases))]
        public void ConstructorTestWithDefaultSettings(MetaMorpheusTask task)
        {
            Directory.CreateDirectory(Path.Combine(GlobalVariables.DataDir, "DefaultParameters"));
            string oldDefaultName = task switch
            {
                SearchTask search => "SearchTaskDefault.toml",
                CalibrationTask calibration => "CalibrationTaskDefault.toml",
                GptmdTask gptmd => "GptmdTaskDefault.toml",
                XLSearchTask xLSearch => "XLSearchTaskDefault.toml",
                GlycoSearchTask glycoSearch => "GlycoSearchTaskDefault.toml",
                SpectralAveragingTask spectralAverage => "SpectralAverageTaskDefault.toml",
                _ => "",
            };
            String pathToCheck = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "DefaultParameters\\" + oldDefaultName);
            Toml.WriteFile<MetaMorpheusTask>(task, pathToCheck, MetaMorpheusTask.tomlConfig);
            var fakeTaskWindow = new FakeGuiWindow(task.TaskType);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(task,
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var testDictionary = testTaskSettingViewModel.AllSettingsDict;
            String defaultName = task switch
            {
                SearchTask search => "SearchTaskDefault(Default)",
                CalibrationTask calibration => "CalibrationTaskDefault(Default)",
                GptmdTask gptmd => "GptmdTaskDefault(Default)",
                XLSearchTask xLSearch => "XLSearchTaskDefault(Default)",
                GlycoSearchTask glycoSearch => "GlycoSearchTaskDefault(Default)",
                SpectralAveragingTask spectralAverage => "SpectralAverageTaskDefault(Default)",
                _ => "",
            };
            Assert.That(testDictionary.Keys.Contains(defaultName));
            var testWindow = fakeTaskWindow.GetTaskFromGui();
            var comparingWindow = testDictionary[defaultName];
            Assert.That(!testWindow.Equals(null));
            TomlFileFolderSerializer.Save("testFile", testWindow);
            string pathOfTestFile = TomlFileFolderSerializer.GetFilePath(testWindow.GetType(), "testFile");
            string pathOfComparingFile =
                TomlFileFolderSerializer.GetFilePath(comparingWindow.GetType(), defaultName);
            var linesA = File.ReadAllLines(pathOfTestFile);
            var linesB = File.ReadAllLines(pathOfComparingFile);

            Assert.That(linesA.Length, Is.EqualTo(linesB.Length));
            for (int i = 0; i < linesA.Length; i++)
            {
                Assert.That(linesA[i], Is.EqualTo(linesB[i]));
            }
        }

        /// <summary>
        /// test the function of saveSettings when modifying non-default task
        /// </summary>
        [Test]
        public void SaveSettingsTest()
        {
            //Directory Set up
            var tomlFile1 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"),
                MetaMorpheusTask.tomlConfig);
            TomlFileFolderSerializer.Save("customBCZ", tomlFile1);
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(),
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            testTaskSettingViewModel.SelectedSettings = "customBCZ";
            testTaskSettingViewModel.SaveSettings();
            var currentTask = testTaskSettingViewModel.AllSettingsDict[testTaskSettingViewModel.SelectedSettings];
            Assert.That(currentTask.CommonParameters.TotalPartitions, Is.EqualTo(1));

            // test change in window is reflected when saving settings
            currentTask = testTaskSettingViewModel.AllSettingsDict[testTaskSettingViewModel.SelectedSettings];
            Assert.That(currentTask.CommonParameters.TotalPartitions, Is.EqualTo(1));
            fakeTaskWindow.TheTask.CommonParameters.TotalPartitions = 36;
            testTaskSettingViewModel.SaveSettings();

            currentTask = testTaskSettingViewModel.AllSettingsDict[testTaskSettingViewModel.SelectedSettings];
            Assert.That(currentTask.CommonParameters.TotalPartitions, Is.EqualTo(36));

            // read back in saved file and check if total partitions was saved
            testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(),
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            testTaskSettingViewModel.SelectedSettings = "customBCZ";
            currentTask = testTaskSettingViewModel.AllSettingsDict[testTaskSettingViewModel.SelectedSettings];
            Assert.That(currentTask.CommonParameters.TotalPartitions, Is.EqualTo(36));
        }

        /// <summary>
        /// test the function of saveSettings when modifying default task
        /// </summary>
        [Test]
        public void SaveDefaultSettingsTest()
        {
            //Directory Set up
            var tomlFile1 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"),
                MetaMorpheusTask.tomlConfig);
            TomlFileFolderSerializer.Save("customBCZ", tomlFile1);
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(),
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);

            //Try to save a default setting
            fakeTaskWindow.TheTask.CommonParameters.TotalPartitions = 36;
            testTaskSettingViewModel.SelectedSettings = "customBCZ";
            testTaskSettingViewModel.SaveAsDefaultSettings();
            testTaskSettingViewModel.SelectedSettings = "customBCZ(Default)";
            testTaskSettingViewModel.SaveSettings();
            testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(),
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var a = testTaskSettingViewModel.AllSettingsDict;
            testTaskSettingViewModel.SelectedSettings = "customBCZ(Default)";
            var currentTask = testTaskSettingViewModel.AllSettingsDict[testTaskSettingViewModel.SelectedSettings];
            Assert.That(currentTask.CommonParameters.TotalPartitions, Is.EqualTo(1));
        }

        /// <summary>
        /// test the function of saveSettingsFromWindow
        /// </summary>
        [Test]
        public void SaveSettingsFromWindowTest()
        {
            //Directory Set up
            var tomlFile1 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"),
                MetaMorpheusTask.tomlConfig);
            TomlFileFolderSerializer.Save("customBCZ", tomlFile1);
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(),
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var testDictionary = testTaskSettingViewModel.AllSettingsDict;

            //Save a task from the window popped up
            testTaskSettingViewModel.SelectedSettings = "customBCZ";
            testTaskSettingViewModel.SaveSettingsFromWindow();//
            Assert.That(!testDictionary.ContainsKey("testTask"));
            testTaskSettingViewModel.TypedSettingsName = "customBCZ";
            testTaskSettingViewModel.SaveSettingsFromWindow();
            Assert.That(!testDictionary.ContainsKey("testTask"));
            testTaskSettingViewModel.TypedSettingsName = "testTask";
            testTaskSettingViewModel.SaveSettingsFromWindow();
            Assert.That(testDictionary.ContainsKey("testTask"));
        }

        /// <summary>
        /// test the function of saveAsDefaultSetting
        /// </summary>
        [Test]
        public void SaveAsDefaultSettingsTest()
        {
            //Directory Set up
            var tomlFile1 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"),
                MetaMorpheusTask.tomlConfig);
            TomlFileFolderSerializer.Save("customBCZ", tomlFile1);
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(),
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var testDictionary = testTaskSettingViewModel.AllSettingsDict;

            //When SelectedSettings == null
            testTaskSettingViewModel.SaveAsDefaultSettings();
            Assert.That(testDictionary.ContainsKey("customBCZ"));
            testTaskSettingViewModel.SelectedSettings = "customBCZ";
            testTaskSettingViewModel.SaveAsDefaultSettings();
            Assert.That(testDictionary.ContainsKey("customBCZ(Default)"));

            //Setting default when there is an existing default
            var tomlFile2 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBY.toml"),
                MetaMorpheusTask.tomlConfig);
            TomlFileFolderSerializer.Save("customBY", tomlFile2); 
            fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(),
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            testDictionary = testTaskSettingViewModel.AllSettingsDict;
            testTaskSettingViewModel.SelectedSettings = "customBY";
            testTaskSettingViewModel.SaveAsDefaultSettings();
            Assert.That(testDictionary.ContainsKey("customBY(Default)"));
            Assert.That(!testDictionary.ContainsKey("customBCZ(Default)"));
        }

        /// <summary>
        /// test the function of deleteSettings when nothing is selected
        /// </summary>
        [Test]
        public void DeleteSettingsTestWithoutSelectedSettings()
        {
            //Directory Set up
            var tomlFile1 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"),
                MetaMorpheusTask.tomlConfig);
            var tomlFile2 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBY.toml"),
                MetaMorpheusTask.tomlConfig);
            var tomlFile3 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customCZ.toml"),
                MetaMorpheusTask.tomlConfig);
            TomlFileFolderSerializer.Save("customBCZ.toml", tomlFile1);
            TomlFileFolderSerializer.Save("customBY.toml", tomlFile2);
            TomlFileFolderSerializer.Save("customCZ.toml", tomlFile3);
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(),
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var testDictionary = testTaskSettingViewModel.AllSettingsDict;

            var expected = new List<string>() { "customBCZ.toml", "customBY.toml", "customCZ.toml" };
            CollectionAssert.AreEqual(expected, testTaskSettingViewModel.AllSettings);
            testTaskSettingViewModel.TypedSettingsName = "Name";
            Assert.That(testTaskSettingViewModel.TypedSettingsName, Is.EqualTo("Name"));

            //test
            Type deletedFileType = tomlFile1.GetType();
            String path = TomlFileFolderSerializer.GetFilePath(deletedFileType, "customBCZ.toml");
            testTaskSettingViewModel.DeleteSettings();
            Assert.That(File.Exists(path), Is.True);
        }

        /// <summary>
        /// test the function of deleteSettings when default setting is selected
        /// </summary>
        [Test]
        public void DeleteSettingsTestSelectingDefault()
        {
            //Directory Set up
            var tomlFile1 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"),
                MetaMorpheusTask.tomlConfig);
            var tomlFile2 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBY.toml"),
                MetaMorpheusTask.tomlConfig);
            var tomlFile3 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customCZ.toml"),
                MetaMorpheusTask.tomlConfig);
            TomlFileFolderSerializer.Save("customBCZ", tomlFile1);
            TomlFileFolderSerializer.Save("customBY", tomlFile2);
            TomlFileFolderSerializer.Save("customCZ", tomlFile3);
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(),
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var testDictionary = testTaskSettingViewModel.AllSettingsDict;


            // save customBCZ as default settings
            testTaskSettingViewModel.SelectedSettings = "customBCZ";
            testTaskSettingViewModel.SaveAsDefaultSettings();
            testTaskSettingViewModel.SelectedSettings = "customBCZ(Default)";

            // try to delete default settings
            testTaskSettingViewModel.DeleteSettings();

            // ensure default was not deleted
            Assert.That(testDictionary.ContainsKey("customBCZ(Default)"));
        }

        /// <summary>
        /// test the function of deleteSettings when select a non-default task
        /// </summary>
        [Test]
        public void DeleteSettingsTestWithSelectedSettingsAndNotDefault()
        {
            //Directory Set up
            var tomlFile1 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"),
                MetaMorpheusTask.tomlConfig);
            var tomlFile2 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBY.toml"),
                MetaMorpheusTask.tomlConfig);
            var tomlFile3 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customCZ.toml"),
                MetaMorpheusTask.tomlConfig);
            TomlFileFolderSerializer.Save("customBCZ.toml", tomlFile1);
            TomlFileFolderSerializer.Save("customBY.toml", tomlFile2);
            TomlFileFolderSerializer.Save("customCZ.toml", tomlFile3);
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(),
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var testDictionary = testTaskSettingViewModel.AllSettingsDict;

            //test
            Type deletedFileType = tomlFile1.GetType();
            String path = TomlFileFolderSerializer.GetFilePath(deletedFileType, "customBCZ.toml");
            testTaskSettingViewModel.SelectedSettings = "customBCZ.toml";
            testTaskSettingViewModel.DeleteSettings();
            Assert.That(!testDictionary.ContainsKey("customBCZ.toml"));
            Assert.That(!File.Exists(path), Is.True);
        }

        [Test]
        public void InstanceTest()
        {
            //Directory Set up
            var tomlFile1 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"),
                MetaMorpheusTask.tomlConfig);
            var tomlFile2 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBY.toml"),
                MetaMorpheusTask.tomlConfig);
            var tomlFile3 = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customCZ.toml"),
                MetaMorpheusTask.tomlConfig);
            TomlFileFolderSerializer.Save("customBCZ.toml", tomlFile1);
            TomlFileFolderSerializer.Save("customBY.toml", tomlFile2);
            TomlFileFolderSerializer.Save("customCZ.toml", tomlFile3);
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(),
                fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var testDictionary = testTaskSettingViewModel.AllSettingsDict;
            //AllSettings
            var expected = new List<string>() { "customBCZ.toml", "customBY.toml", "customCZ.toml" };
            CollectionAssert.AreEqual(expected, testTaskSettingViewModel.AllSettings);
            //TypedSettingsName
            testTaskSettingViewModel.TypedSettingsName = "Name";
            Assert.That(testTaskSettingViewModel.TypedSettingsName, Is.EqualTo("Name"));
            //SelectedSettings
            testTaskSettingViewModel.SelectedSettings = "customBCZ.toml";
            Assert.That(testTaskSettingViewModel.SelectedSettings.Equals("customBCZ.toml"));
            //AllSettingsDict
            Assert.That(testTaskSettingViewModel.AllSettingsDict.Equals(testDictionary));
            //TheTask
            Assert.That(testTaskSettingViewModel.TheTask.Equals(testDictionary["customBCZ.toml"]));
        }
    }
}
