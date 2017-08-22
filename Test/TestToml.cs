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

            Assert.AreEqual(searchTask.DeconvolutionMassTolerance.ToString(), searchTaskLoaded.DeconvolutionMassTolerance.ToString());
            Assert.AreEqual(searchTask.ProductMassTolerance.ToString(), searchTaskLoaded.ProductMassTolerance.ToString());
            Assert.AreEqual(searchTask.MassDiffAcceptors[0].FileNameAddition, searchTaskLoaded.MassDiffAcceptors[0].FileNameAddition);
            Assert.AreEqual(searchTask.ListOfModsFixed[0].Item1, searchTaskLoaded.ListOfModsFixed[0].Item1);
            Assert.AreEqual(searchTask.ListOfModsFixed[0].Item2, searchTaskLoaded.ListOfModsFixed[0].Item2);
            Assert.AreEqual(searchTask.ListOfModsLocalize.Count, searchTaskLoaded.ListOfModsLocalize.Count);
            Assert.AreEqual(searchTask.ListOfModsFixed.Count, searchTaskLoaded.ListOfModsFixed.Count);
            Assert.AreEqual(searchTask.ListOfModsVariable.Count, searchTaskLoaded.ListOfModsVariable.Count);

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
            SearchTask searchTask = new SearchTask();
            var dir = Directory.GetParent(Directory.GetParent(Directory.GetCurrentDirectory()).ToString()).ToString();
            Console.WriteLine();
            var file = Directory.GetFiles(dir.ToString(), Path.GetFileNameWithoutExtension("testFileSpecfic") + ".to*");
            var fileSpecificToml = Toml.ReadFile(file[0], MetaMorpheusTask.tomlConfig);
            var tomlSettingsList = fileSpecificToml.ToDictionary(p => p.Key);
            Assert.AreEqual(tomlSettingsList["Protease"].Value.Get<string>(), "TestCustomProtease");

        }

        #endregion Public Methods
    }


}