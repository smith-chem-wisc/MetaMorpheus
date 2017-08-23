using Nett;
using NUnit.Framework;
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

            Assert.AreEqual(searchTask.commonParameters.DeconvolutionMassTolerance.ToString(), searchTaskLoaded.commonParameters.DeconvolutionMassTolerance.ToString());
            Assert.AreEqual(searchTask.commonParameters.ProductMassTolerance.ToString(), searchTaskLoaded.commonParameters.ProductMassTolerance.ToString());
            Assert.AreEqual(searchTask.searchParameters.MassDiffAcceptors[0].FileNameAddition, searchTaskLoaded.searchParameters.MassDiffAcceptors[0].FileNameAddition);
            Assert.AreEqual(searchTask.commonParameters.ListOfModsFixed[0].Item1, searchTaskLoaded.commonParameters.ListOfModsFixed[0].Item1);
            Assert.AreEqual(searchTask.commonParameters.ListOfModsFixed[0].Item2, searchTaskLoaded.commonParameters.ListOfModsFixed[0].Item2);
            Assert.AreEqual(searchTask.commonParameters.ListOfModsLocalize.Count, searchTaskLoaded.commonParameters.ListOfModsLocalize.Count);
            Assert.AreEqual(searchTask.commonParameters.ListOfModsFixed.Count, searchTaskLoaded.commonParameters.ListOfModsFixed.Count);
            Assert.AreEqual(searchTask.commonParameters.ListOfModsVariable.Count, searchTaskLoaded.commonParameters.ListOfModsVariable.Count);

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

        #endregion Public Methods
    }
}