using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Security.AccessControl;
using System.Text;
using System.Threading.Tasks;
using EngineLayer;
using GuiFunctions;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using TaskLayer;

namespace Test
{
    public class TaskSettingViewModelTests
    {
        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, "SettingsDataFiles"));
            var tomlFile1 = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\customBCZ.toml"), MetaMorpheusTask.tomlConfig);
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

        [Test]
        public void testSaveSettings()
        {
            TaskSettingViewModel testTaskSettingViewModel = new TaskSettingViewModel(new SearchTask());

            var paramsOne = new CommonParameters(null, DissociationType.CID);
            var paramsTwo = new CommonParameters(null, DissociationType.HCD);
            var paramsThree = new CommonParameters(null, DissociationType.PD);

            var one = new SearchTask();
            var two = new SearchTask();
            var three = new SearchTask();

            one.CommonParameters = paramsOne;
            two.CommonParameters = paramsTwo;
            three.CommonParameters = paramsThree;

        }
    }
}
