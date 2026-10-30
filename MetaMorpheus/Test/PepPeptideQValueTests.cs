using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// PEP trains at the peptide level once there are 100 distinct full sequences, and labels its positives by
    /// peptide q-value. When there are more PSMs than <see cref="FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue"/>
    /// but not more peptides, those peptide q-values must still be computed before PEP, or no target is a positive.
    /// </summary>
    [TestFixture]
    [NonParallelizable] // sets FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue, a process-wide static
    public static class PepPeptideQValueTests
    {
        [Test]
        public static void PepTrainsWhenThereAreMorePsmsThanTheThresholdButFewerPeptides()
        {
            string testDirectory = TestContext.CurrentContext.TestDirectory;
            string outputFolder = Path.Combine(testDirectory, "PepPeptideQValueTests");
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);

            var task = Nett.Toml.ReadFile<SearchTask>(Path.Combine(testDirectory, @"TestData\Task1-SearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            task.CommonParameters.QValueCutoffForPepCalculation = 0.01;
            var spectraFiles = new List<string>
            {
                Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip.mzML"),
                Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip_2.mzML")
            };
            var databases = new List<DbForTask> { new DbForTask(Path.Combine(testDirectory, @"TestData\TaGe_SA_A549_3_snip.fasta"), false) };

            var threshold = typeof(FdrAnalysisEngine).GetProperty(nameof(FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue))!;
            int original = FdrAnalysisEngine.PsmCountThresholdForInvertedQvalue;
            lock (EverythingRunnerEngineTestCase.myLock)
            {
                try
                {
                    threshold.SetValue(null, 500); // 770 PSMs, 308 peptides
                    new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("search", task) }, spectraFiles, databases, outputFolder).Run();
                }
                finally
                {
                    threshold.SetValue(null, original);
                }
            }

            string results = File.ReadAllText(Path.Combine(outputFolder, "search", "results.txt"));
            Assert.That(results, Does.Not.Contain("Posterior error probability analysis failed"));
            Assert.That(results, Does.Contain("Metrics for Determination of PEP"));
            Directory.Delete(outputFolder, true);
        }
    }
}
