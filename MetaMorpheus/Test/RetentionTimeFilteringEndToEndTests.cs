using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;

namespace Test
{
    [TestFixture]
    [NonParallelizable]
    public class RetentionTimeFilteringEndToEndTests
    {
        private static readonly FieldInfo ClassicSearchMatches = typeof(ClassicSearchEngine)
            .GetField("SpectralMatches", BindingFlags.Instance | BindingFlags.NonPublic);

        [Test]
        public void SearchTask_RestrictsReportedPsmRetentionTimes()
        {
            string testDirectory = TestContext.CurrentContext.TestDirectory;
            string configPath = Path.Combine(testDirectory, "TestData", "Task1-SearchTaskconfig.toml");
            string spectraPath = Path.Combine(testDirectory, "TestData", "TaGe_SA_A549_3_snip.mzML");
            string databasePath = Path.Combine(testDirectory, "TestData", "TaGe_SA_A549_3_snip.fasta");

            List<double> unrestricted = RunSearchTask<SearchTask>(configPath, spectraPath, databasePath, null, null);
            AssertRetentionRangeCanBeTested(unrestricted);
            (double minimum, double maximum) = SelectInteriorRetentionRange(unrestricted);

            List<double> restricted = RunSearchTask<SearchTask>(configPath, spectraPath, databasePath, minimum, maximum);

            AssertRestrictedResults(unrestricted, restricted, minimum, maximum);
        }

        [Test]
        public void GptmdTask_RestrictsClassicSearchRetentionTimes()
        {
            string testDirectory = TestContext.CurrentContext.TestDirectory;
            string configPath = Path.Combine(testDirectory, "TestData", "ModificationAnalysis", "GPTMDTaskconfig.toml");
            string spectraPath = Path.Combine(testDirectory, "TestData", "ModificationAnalysis", "modificationAnalysis.mzML");
            string databasePath = Path.Combine(testDirectory, "TestData", "ModificationAnalysis", "modificationAnalysis.fasta");

            List<double> unrestricted = RunSearchTask<GptmdTask>(configPath, spectraPath, databasePath, null, null);
            AssertRetentionRangeCanBeTested(unrestricted);
            (double minimum, double maximum) = SelectInteriorRetentionRange(unrestricted);

            List<double> restricted = RunSearchTask<GptmdTask>(configPath, spectraPath, databasePath, minimum, maximum);

            AssertRestrictedResults(unrestricted, restricted, minimum, maximum);
        }

        private static List<double> RunSearchTask<TTask>(string configPath, string spectraPath, string databasePath,
            double? minimumRetentionTime, double? maximumRetentionTime) where TTask : MetaMorpheusTask
        {
            string outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "RetentionTimeFiltering", typeof(TTask).Name + "_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(outputDirectory);

            string taskConfigPath = CreateTaskConfig(configPath, outputDirectory, minimumRetentionTime, maximumRetentionTime);
            TTask task = Toml.ReadFile<TTask>(taskConfigPath, MetaMorpheusTask.tomlConfig);
            List<double> retentionTimes = new();

            EventHandler<SingleEngineFinishedEventArgs> finishedEngineHandler = (sender, args) =>
            {
                if (args.MyResults.MyEngine is not ClassicSearchEngine classicSearchEngine)
                    return;

                SpectralMatch[] matches = (SpectralMatch[])ClassicSearchMatches.GetValue(classicSearchEngine);
                retentionTimes.AddRange(matches.Where(match => match != null).Select(match => match.ScanRetentionTime));
            };

            MetaMorpheusEngine.FinishedSingleEngineHandler += finishedEngineHandler;
            try
            {
                task.RunTask(outputDirectory,
                    new List<DbForTask> { new DbForTask(databasePath, false) },
                    new List<string> { spectraPath },
                    typeof(TTask).Name);

                if (typeof(TTask) == typeof(SearchTask))
                {
                    string psmPath = Path.Combine(outputDirectory, "AllPSMs.psmtsv");
                    var psms = SpectrumMatchTsvReader.ReadPsmTsv(psmPath, out List<string> warnings);
                    retentionTimes = psms.Where(psm => psm != null).Select(psm => psm.RetentionTime).ToList();
                }
            }
            finally
            {
                MetaMorpheusEngine.FinishedSingleEngineHandler -= finishedEngineHandler;
                if (File.Exists(taskConfigPath))
                    File.Delete(taskConfigPath);
                if (Directory.Exists(outputDirectory))
                    Directory.Delete(outputDirectory, true);
            }

            return retentionTimes;
        }

        private static string CreateTaskConfig(string sourceConfigPath, string outputDirectory,
            double? minimumRetentionTime, double? maximumRetentionTime)
        {
            string taskConfigPath = Path.Combine(outputDirectory, "task.toml");
            List<string> lines = File.ReadAllLines(sourceConfigPath).ToList();

            if (minimumRetentionTime.HasValue || maximumRetentionTime.HasValue)
            {
                int digestionSection = lines.FindIndex(line => line.Equals("[CommonParameters.DigestionParams]", StringComparison.Ordinal));
                Assert.That(digestionSection, Is.GreaterThan(0), "The test task must contain CommonParameters.DigestionParams.");
                string maximum = maximumRetentionTime.HasValue
                    ? maximumRetentionTime.Value.ToString(System.Globalization.CultureInfo.InvariantCulture)
                    : "MaxValue";
                lines.Insert(digestionSection, $"RetentionTimeRange = \"{(minimumRetentionTime ?? 0).ToString(System.Globalization.CultureInfo.InvariantCulture)};{maximum}\"");
            }

            File.WriteAllLines(taskConfigPath, lines);
            return taskConfigPath;
        }

        private static void AssertRetentionRangeCanBeTested(List<double> retentionTimes)
        {
            Assert.That(retentionTimes, Is.Not.Empty);
            Assert.That(retentionTimes.Distinct().Count(), Is.GreaterThanOrEqualTo(4),
                "The fixture must contain enough distinct PSM retention times to test a restricted interval.");
        }

        private static (double Minimum, double Maximum) SelectInteriorRetentionRange(List<double> retentionTimes)
        {
            List<double> ordered = retentionTimes.Distinct().OrderBy(retentionTime => retentionTime).ToList();
            return (ordered[ordered.Count / 4], ordered[(ordered.Count * 3) / 4]);
        }

        private static void AssertRestrictedResults(List<double> unrestricted, List<double> restricted,
            double minimumRetentionTime, double maximumRetentionTime)
        {
            Assert.That(restricted, Is.Not.Empty);
            Assert.That(restricted.Count, Is.LessThan(unrestricted.Count));
            Assert.That(unrestricted.Any(retentionTime => retentionTime < minimumRetentionTime || retentionTime > maximumRetentionTime), Is.True);
            Assert.That(restricted, Is.All.GreaterThanOrEqualTo(minimumRetentionTime));
            Assert.That(restricted, Is.All.LessThanOrEqualTo(maximumRetentionTime));
        }
    }
}
