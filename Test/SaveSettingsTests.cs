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
            TomlFileFolderSerializer.Save("customBCZ.toml", tomlFile1);
            TomlFileFolderSerializer.Save("customBY.toml", tomlFile2);
            TomlFileFolderSerializer.Save("customCZ.toml", tomlFile3);
        }

        [OneTimeTearDown]
        public static void OneTimeTearDown()
        {
            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "SettingsDataFiles"), true);
        }

        private MetaMorpheusTask GetFakeSearchTaskFromGui() => new SearchTask();
        private MetaMorpheusTask GetFakeGptmdTaskFromGui() => new GptmdTask();
        private MetaMorpheusTask GetFakeCalibrationTaskFromGui() => new CalibrationTask();
        private MetaMorpheusTask GetFakeXLSearchTaskFromGui() => new XLSearchTask();
        private MetaMorpheusTask GetFakeGlycoSearchTaskFromGui() => new GlycoSearchTask();

        private void UpdateFieldsFromNewTask(MetaMorpheusTask task)
        {

        }
        
        [Test]
        public void testConstructor()
        {
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask(), UpdateFieldsFromNewTask, GetFakeSearchTaskFromGui);
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
    }
}
