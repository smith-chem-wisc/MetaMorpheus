using EngineLayer;
using Nett;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public class TestToml
    {
        #region Public Methods

        [Test]
        public static void TestTomlFunction()
        {
            SearchTask searchTask = new SearchTask();
            Toml.WriteFile(searchTask, "SearchTask.toml", MetaMorpheusTask.tomlConfig);
            var searchTaskLoaded = Toml.ReadFile<SearchTask>("SearchTask.toml", MetaMorpheusTask.tomlConfig);

            Assert.AreEqual(searchTask.CommonParameters.DeconvolutionMassTolerance.ToString(), searchTaskLoaded.CommonParameters.DeconvolutionMassTolerance.ToString());
            Assert.AreEqual(searchTask.CommonParameters.ProductMassTolerance.ToString(), searchTaskLoaded.CommonParameters.ProductMassTolerance.ToString());
            Assert.AreEqual(searchTask.SearchParameters.MassDiffAcceptors[0].FileNameAddition, searchTaskLoaded.SearchParameters.MassDiffAcceptors[0].FileNameAddition);
            Assert.AreEqual(searchTask.CommonParameters.ListOfModsFixed[0].Item1, searchTaskLoaded.CommonParameters.ListOfModsFixed[0].Item1);
            Assert.AreEqual(searchTask.CommonParameters.ListOfModsFixed[0].Item2, searchTaskLoaded.CommonParameters.ListOfModsFixed[0].Item2);
            Assert.AreEqual(searchTask.CommonParameters.ListOfModsLocalize, searchTaskLoaded.CommonParameters.ListOfModsLocalize);
            Assert.AreEqual(searchTask.CommonParameters.ListOfModsFixed.Count, searchTaskLoaded.CommonParameters.ListOfModsFixed.Count);
            Assert.AreEqual(searchTask.CommonParameters.ListOfModsVariable.Count, searchTaskLoaded.CommonParameters.ListOfModsVariable.Count);

            CalibrationTask calibrationTask = new CalibrationTask();
            Toml.WriteFile(calibrationTask, "CalibrationTask.toml", MetaMorpheusTask.tomlConfig);
            var calibrationTaskLoaded = Toml.ReadFile<CalibrationTask>("CalibrationTask.toml", MetaMorpheusTask.tomlConfig);

            GptmdTask gptmdTask = new GptmdTask();
            Toml.WriteFile(gptmdTask, "GptmdTask.toml", MetaMorpheusTask.tomlConfig);
            var gptmdTaskLoaded = Toml.ReadFile<GptmdTask>("GptmdTask.toml", MetaMorpheusTask.tomlConfig);

            XLSearchTask xLSearchTask = new XLSearchTask();
            Toml.WriteFile(xLSearchTask, "XLSearchTask.toml", MetaMorpheusTask.tomlConfig);
            var xLSearchTaskLoaded = Toml.ReadFile<XLSearchTask>("XLSearchTask.toml", MetaMorpheusTask.tomlConfig);
        }

        [Test]
        public static void TestTomlForSpecficFiles()
        {
            var dir = Directory.GetParent(Directory.GetParent(Directory.GetCurrentDirectory()).ToString()).ToString();
            var file = Directory.GetFiles(dir.ToString(), Path.GetFileNameWithoutExtension("testFileSpecfic") + ".to*");
            var fileSpecificToml = Toml.ReadFile(file[0], MetaMorpheusTask.tomlConfig);
            var tomlSettingsList = fileSpecificToml.ToDictionary(p => p.Key);
            Assert.AreEqual(tomlSettingsList["Protease"].Value.Get<string>(), "TestCustomProtease");

        }

        #endregion Public Methods
    }


}