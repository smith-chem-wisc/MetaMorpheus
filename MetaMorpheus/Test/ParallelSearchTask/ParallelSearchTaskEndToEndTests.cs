using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer.DatabaseLoading;
using Nett;
using NUnit.Framework;
using TaskLayer;
using TaskLayer.ParallelSearch;

namespace Test.ParallelSearchTask;

[TestFixture]
public class ParallelSearchTaskEndToEndTests
{
    private static readonly string TestDataFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
    private static readonly string SpectraPath = Path.Combine(TestDataFolder, "TaGe_SA_A549_3_snip.mzML");
    private static readonly string BaseDatabasePath = Path.Combine(TestDataFolder, "TaGe_SA_A549_3_snip.fasta");
    private static readonly string[] TransientDatabasePaths =
    [
        Path.Combine(TestDataFolder, "hela_snip_for_unitTest.fasta"),
        Path.Combine(TestDataFolder, "smalldb.fasta")
    ];

    private static TaskLayer.ParallelSearch.ParallelSearchTask CreateTask(bool overwrite)
    {
        var task = new TaskLayer.ParallelSearch.ParallelSearchTask(
            TransientDatabasePaths.Select(p => new DbForTask(p, false)).ToList());
        task.ParallelSearchParameters.OverwriteTransientSearchOutputs = overwrite;
        task.ParallelSearchParameters.MaxSearchesInParallel = 2;
        return task;
    }

    private static Dictionary<string, TaskLayer.ParallelSearch.Analysis.TransientDatabaseMetrics> ReadSummary(string outputFolder)
    {
        var cache = new ParallelSearchResultCache(
            Path.Combine(outputFolder, TransientDatabaseResultsManager.SummaryResultsFileName));
        cache.InitializeCache();
        return cache.AllResultsDictionary;
    }

    [Test]
    public void RunTask_WritesOneSummaryRowPerTransientDatabase_AndReusesThemOnRerun()
    {
        string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "ParallelSearchEndToEnd");
        if (Directory.Exists(outputFolder))
            Directory.Delete(outputFolder, true);
        Directory.CreateDirectory(outputFolder);

        try
        {
            var baseDatabases = new List<DbForTask> { new(BaseDatabasePath, false) };
            var spectra = new List<string> { SpectraPath };
            string[] expectedNames = TransientDatabasePaths.Select(Path.GetFileNameWithoutExtension).ToArray();

            CreateTask(overwrite: true).RunTask(outputFolder, baseDatabases, spectra, "parallel");

            string summaryPath = Path.Combine(outputFolder, TransientDatabaseResultsManager.SummaryResultsFileName);
            Assert.That(File.Exists(summaryPath), Is.True);

            string[] summaryLines = File.ReadAllLines(summaryPath);
            string[] header = summaryLines[0].Split(',');
            var firstRun = ReadSummary(outputFolder);

            Assert.Multiple(() =>
            {
                // header plus one row per transient database
                Assert.That(summaryLines.Length, Is.EqualTo(expectedNames.Length + 1));
                Assert.That(header, Does.Contain("DatabaseName"));
                Assert.That(header, Does.Contain("TransientProteinCount"));
                Assert.That(header, Does.Contain("TargetPsmsFromTransientDbAtQValueThreshold"));
                Assert.That(firstRun.Keys, Is.EquivalentTo(expectedNames));
                Assert.That(firstRun.Values.All(m => m.TransientProteinCount > 0), Is.True);
                Assert.That(firstRun.Values.All(m => m.TotalProteins > m.TransientProteinCount), Is.True);
                Assert.That(File.Exists(Path.Combine(outputFolder, "BasePSMs.psmtsv")), Is.True);
                foreach (string name in expectedNames)
                {
                    Assert.That(File.Exists(Path.Combine(outputFolder, name, $"{name}_AllPSMs.psmtsv")), Is.True);
                    Assert.That(File.Exists(Path.Combine(outputFolder, name, "results.txt")), Is.True);
                }
            });

            // Rerun without overwriting: every database is already in the summary, so nothing is searched again
            string perDatabaseFile = Path.Combine(outputFolder, expectedNames[0], $"{expectedNames[0]}_AllPSMs.psmtsv");
            var writeTimeBeforeRerun = File.GetLastWriteTimeUtc(perDatabaseFile);

            CreateTask(overwrite: false).RunTask(outputFolder, baseDatabases, spectra, "parallel");

            var secondRun = ReadSummary(outputFolder);
            Assert.Multiple(() =>
            {
                Assert.That(File.GetLastWriteTimeUtc(perDatabaseFile), Is.EqualTo(writeTimeBeforeRerun));
                Assert.That(File.ReadAllLines(summaryPath).Length, Is.EqualTo(expectedNames.Length + 1));
                Assert.That(secondRun.Keys, Is.EquivalentTo(expectedNames));
                foreach (string name in expectedNames)
                {
                    Assert.That(secondRun[name].TransientProteinCount, Is.EqualTo(firstRun[name].TransientProteinCount));
                    Assert.That(secondRun[name].TargetPsmsFromTransientDbAtQValueThreshold,
                        Is.EqualTo(firstRun[name].TargetPsmsFromTransientDbAtQValueThreshold));
                }
            });
        }
        finally
        {
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);
        }
    }

    [Test]
    public void Toml_RoundTripsTaskTypeAndTransientDatabases()
    {
        string tomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "ParallelSearchRoundTrip.toml");
        try
        {
            var task = CreateTask(overwrite: false);
            Toml.WriteFile(task, tomlPath, MetaMorpheusTask.tomlConfig);

            // the command line picks the task class from this value
            string taskType = Toml.ReadFile(tomlPath, MetaMorpheusTask.tomlConfig).Get<string>("TaskType");
            var readBack = Toml.ReadFile<TaskLayer.ParallelSearch.ParallelSearchTask>(tomlPath, MetaMorpheusTask.tomlConfig);

            Assert.Multiple(() =>
            {
                Assert.That(taskType, Is.EqualTo(nameof(MyTask.ParallelSearch)));
                Assert.That(readBack.TaskType, Is.EqualTo(MyTask.ParallelSearch));
                Assert.That(readBack.ParallelSearchParameters.TransientDatabases.Select(db => db.FilePath),
                    Is.EqualTo(TransientDatabasePaths));
                Assert.That(readBack.ParallelSearchParameters.OverwriteTransientSearchOutputs, Is.False);
                Assert.That(readBack.ParallelSearchParameters.MaxSearchesInParallel, Is.EqualTo(2));
                Assert.That(readBack.SearchParameters, Is.SameAs(readBack.ParallelSearchParameters));
            });
        }
        finally
        {
            if (File.Exists(tomlPath))
                File.Delete(tomlPath);
        }
    }
}
