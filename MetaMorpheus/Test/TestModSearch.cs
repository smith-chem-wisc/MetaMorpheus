using Easy.Common.Extensions;
using Nett;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using EngineLayer.GlycoSearch;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public class TestModSearch
    {
        [Test]
        public static void TestModSearchEngine()
        {
            string outputFolder = "";

            var glycoSearchTask = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\GlycoSearchTaskconfigOGlycoTest_Run.toml"), MetaMorpheusTask.tomlConfig);
            glycoSearchTask._glycoSearchParameters.GlycoSearchType = GlycoSearchType.ModSearch;
            DbForTask db = new(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\P16150.fasta"), false);
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"GlycoTestData\2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf");
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", glycoSearchTask) }, new List<string> { spectraFile }, new List<DbForTask> { db }, outputFolder).Run();
        }
    }
}
