using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Security.AccessControl;
using System.Text;
using System.Threading.Tasks;
using GuiFunctions;
using Nett;
using NUnit.Framework;
using TaskLayer;

namespace Test
{

    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SaveSettingsTests
    {
        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, "SettingsDataFiles"));
            var tomlFile1 = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory,@"TestData\customBCZ.toml"), MetaMorpheusTask.tomlConfig);
            var tomlFile2 = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBY.toml"), MetaMorpheusTask.tomlConfig);
            var tomlFile3 = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customCZ.toml"), MetaMorpheusTask.tomlConfig);
            TomlFileFolderSerializer.Save("customBCZ", tomlFile1);
            TomlFileFolderSerializer.Save("customBY", tomlFile2);
            TomlFileFolderSerializer.Save("customCZ", tomlFile3);

            // TODO: Repeat above for callibration task, without setting a default

            // example code
            var task = new CalibrationTask();
            Toml.WriteFile<CalibrationTask>(task, "file out path.toml", MetaMorpheusTask.tomlConfig);
        }   

        [OneTimeTearDown]
        public static void OneTimeTearDown()
        {
            //Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "SettingsDataFiles"), true);
        }

        private class FakeGuiWindow
        {
            public MetaMorpheusTask TheTask { get; private set; }
            public FakeGuiWindow(MyTask taskType)
            {
                TheTask = taskType switch
                {
                    MyTask.Search => new SearchTask(),
                    MyTask.Gptmd => new GptmdTask(),
                    MyTask.Calibrate => new CalibrationTask(),
                    MyTask.XLSearch => new XLSearchTask(),
                    MyTask.GlycoSearch => new GlycoSearchTask(),
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

        // TODO: Test what happens when no default is set

        [Test]
        public void ConstructorTestWithoutDefaultSettings()
        {
            // set up directories



            // test
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(), fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var testDictionary = testTaskSettingViewModel.AllSettingsDict;
            Assert.That(testDictionary.ContainsKey("DefaultSetting(Default)"));
            var testWindow = fakeTaskWindow.GetTaskFromGui();
            var comparingWindow = testDictionary["DefaultSetting(Default)"];
            Assert.That(testWindow.Equals(comparingWindow));
            TomlFileFolderSerializer.Save("testFile", testWindow);
            string pathOfTestFile = TomlFileFolderSerializer.GetFilePath(testWindow.GetType(), "testFile");
            string pathOfComparingFile = TomlFileFolderSerializer.GetFilePath(comparingWindow.GetType(), "DefaultSetting(Default)");
            var linesA = File.ReadAllLines(pathOfTestFile);
            var linesB = File.ReadAllLines(pathOfComparingFile);

            Assert.That(linesA.Length, Is.EqualTo(linesB.Length));
            for (int i = 0; i < linesA.Length; i++)
            {
                //if (linesA[i].StartsWith("MaxMissedCleavages"))
                //{
                //    Assert.That(linesA[i], Is.Not.EqualTo(linesB[i]));
                //    Assert.That(linesA[i].EndsWith('1'));
                //    Assert.That(linesB[i].EndsWith('2'));
                //}
                //else
                Assert.That(linesA[i], Is.EqualTo(linesB[i]));
            }


            // clean up directories
        }

        [Test]
        public void ExampleTest()
        {
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(), fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);

            // test everything is loaded correctly
            // check task fake window is correct
            // check serialized files is correct
            // check dictionary is correct

            // run a command


            // test everything has changed as it should
            // check task fake window is correct
            // check serialized files is correct
            // check dictionary is correct


            // example read all lines comparison
            string pathA = "";
            string pathB = "";

            var linesA = File.ReadAllLines(pathA);
            var linesB = File.ReadAllLines(pathB);

            Assert.That(linesA.Length, Is.EqualTo(linesB.Length));
            for (int i = 0; i < linesA.Length; i++)
            {
                if (linesA[i].StartsWith("MaxMissedCleavages"))
                {
                    Assert.That(linesA[i], Is.Not.EqualTo(linesB[i]));
                    Assert.That(linesA[i].EndsWith('1'));
                    Assert.That(linesB[i].EndsWith('2'));
                }
                else
                    Assert.That(linesA[i], Is.EqualTo(linesB[i]));
            }
            

        }

        [Test]
        public void testConstructor()
        {
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(), fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var a = testTaskSettingViewModel.AllSettings;
            string[] testArrayOfKeys = testTaskSettingViewModel.AllSettingsDict.Keys.ToArray();
            Assert.That(testArrayOfKeys[0].Equals("customBCZ.toml"));
            Assert.That(testArrayOfKeys[1].Equals("customBY.toml"));
            Assert.That(testArrayOfKeys[2].Equals("customCZ.toml"));
            Assert.That(false);
        }


        [Test]
        public void testDeleteGetFilePath()
        {
            var tomlFile1 = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"), MetaMorpheusTask.tomlConfig);
            var tomlFile2 = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBY.toml"), MetaMorpheusTask.tomlConfig);
            var tomlFile3 = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customCZ.toml"), MetaMorpheusTask.tomlConfig);
            TomlFileFolderSerializer.Save("customBCZ.toml", tomlFile1);
            TomlFileFolderSerializer.Save("customBY.toml", tomlFile2);
            TomlFileFolderSerializer.Save("customCZ.toml", tomlFile3);
            Type deletedFileType = tomlFile1.GetType();
            String path = TomlFileFolderSerializer.GetFilePath(deletedFileType, "customBCZ.toml");
            TomlFileFolderSerializer.Delete(deletedFileType, "customBCZ.toml");
            Assert.That(!File.Exists(path),Is.True);                        
        }

        [Test]
        public void testLoadAllNamesOfType()
        {
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask());
            var tomlFile1 = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"), MetaMorpheusTask.tomlConfig);
            Type testFolder = tomlFile1.GetType();
            String[] testStringArray = TomlFileFolderSerializer.LoadAllNamesOfType(testFolder);
            Assert.That(testStringArray[0].Equals("customBCZ.toml"));
            Assert.That(testStringArray[1].Equals("customBY.toml"));
            Assert.That(testStringArray[2].Equals("customCZ.toml"));
        }

        [Test]
        public void testGetFilePath()
        {
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask());
            var tomlFile1 = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"), MetaMorpheusTask.tomlConfig);
            string testPath = TomlFileFolderSerializer.GetFilePath(tomlFile1.GetType(), "customBCZ.toml");
            Assert.That(testPath.Contains("toml"));
        }

        [Test]
        public void testSaveAs()
        {
            var fakeTaskWindow = new FakeGuiWindow(MyTask.Search);
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(), fakeTaskWindow.UpdateFieldsFromNewTask, fakeTaskWindow.GetTaskFromGui);
            var testDictionary = testTaskSettingViewModel.AllSettingsDict;

            // TODO: Replace new search task with refernce B
            var referenceA = TomlFileFolderSerializer.Deserialize<MetaMorpheusTask>("C:/Users/Administrator/source/repos/MetaMorpheus/Test/" +
                "bin/Debug/net6.0-windows/TestData/SavedSettingsReferenceDirectory/a.toml");
            testTaskSettingViewModel.UpdateFieldsInGuiWithNewTask.Invoke(referenceA as MetaMorpheusTask);
            testTaskSettingViewModel.TypedSettingsName = "a";
            testTaskSettingViewModel.SaveSettingsFromWindow();
            Assert.That(false);
        }
    }
}
