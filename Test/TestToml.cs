﻿using EngineLayer;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class TestToml
    {
        #region Public Methods

        [Test]
        public static void TestTomlFunction()
        {
            SearchTask searchTask = new SearchTask
            {
                CommonParameters = new CommonParameters
                {
                    ProductMassTolerance = new PpmTolerance(666),
                    ListOfModsFixed = new List<(string, string)> { ("a", "b"), ("c", "d") },
                    ListOfModsVariable = new List<(string, string)> { ("e", "f"), ("g", "h") },
                    ListOfModTypesLocalize = new List<string> { "i", "j" },
                    LocalizeAll = false,
                },
            };
            Toml.WriteFile(searchTask, "SearchTask.toml", MetaMorpheusTask.tomlConfig);
            var searchTaskLoaded = Toml.ReadFile<SearchTask>("SearchTask.toml", MetaMorpheusTask.tomlConfig);

            Assert.AreEqual(searchTask.CommonParameters.DeconvolutionMassTolerance.ToString(), searchTaskLoaded.CommonParameters.DeconvolutionMassTolerance.ToString());
            Assert.AreEqual(searchTask.CommonParameters.ProductMassTolerance.ToString(), searchTaskLoaded.CommonParameters.ProductMassTolerance.ToString());
            Assert.AreEqual(searchTask.CommonParameters.PrecursorMassTolerance.ToString(), searchTaskLoaded.CommonParameters.PrecursorMassTolerance.ToString());

            Assert.AreEqual(searchTask.CommonParameters.ListOfModsFixed.Count(), searchTaskLoaded.CommonParameters.ListOfModsFixed.Count());
            Assert.AreEqual(searchTask.CommonParameters.ListOfModsFixed.First().Item1, searchTaskLoaded.CommonParameters.ListOfModsFixed.First().Item1);
            Assert.AreEqual(searchTask.CommonParameters.ListOfModsFixed.First().Item2, searchTaskLoaded.CommonParameters.ListOfModsFixed.First().Item2);

            Assert.AreEqual(searchTask.CommonParameters.ListOfModsVariable.Count(), searchTaskLoaded.CommonParameters.ListOfModsVariable.Count());

            Assert.AreEqual(searchTask.CommonParameters.ListOfModTypesLocalize.Count(), searchTaskLoaded.CommonParameters.ListOfModTypesLocalize.Count());
            Assert.AreEqual(searchTask.CommonParameters.ListOfModTypesLocalize.First(), searchTaskLoaded.CommonParameters.ListOfModTypesLocalize.First());
            Assert.AreEqual(searchTask.CommonParameters.ListOfModTypesLocalize.Last(), searchTaskLoaded.CommonParameters.ListOfModTypesLocalize.Last());

            Assert.AreEqual(searchTask.SearchParameters.MassDiffAcceptorType, searchTaskLoaded.SearchParameters.MassDiffAcceptorType);
            Assert.AreEqual(searchTask.SearchParameters.CustomMdac, searchTaskLoaded.SearchParameters.CustomMdac);

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
            var fileSpecificToml = Toml.ReadFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "testFileSpecfic.toml"), MetaMorpheusTask.tomlConfig);
            var tomlSettingsList = fileSpecificToml.ToDictionary(p => p.Key);
            Assert.AreEqual(tomlSettingsList["Protease"].Value.Get<string>(), "Asp-N");
            Assert.IsFalse(tomlSettingsList.ContainsKey("MaxMissedCleavages"));
            Assert.IsFalse(tomlSettingsList.ContainsKey("InitiatorMethionineBehavior"));

            FileSpecificParameters f = new FileSpecificParameters(fileSpecificToml);

            Assert.AreEqual("Asp-N", f.Protease.Name);
            Assert.IsNull(f.MaxMissedCleavages);

            CommonParameters c = MetaMorpheusTask.SetAllFileSpecificCommonParams(new CommonParameters(), f);

            Assert.AreEqual("Asp-N", c.DigestionParams.Protease.Name);
            Assert.AreEqual(2, c.DigestionParams.MaxMissedCleavages);
        }

        #endregion Public Methods
    }
}