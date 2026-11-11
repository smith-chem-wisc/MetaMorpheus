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

    private static string FreshOutputFolder(string name)
    {
        string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
        if (Directory.Exists(outputFolder))
            Directory.Delete(outputFolder, true);
        Directory.CreateDirectory(outputFolder);
        return outputFolder;
    }

    [Test]
    public void RunTask_WithEveryOutputEnabled_WritesAllAndTransientFilesAndCompressesThem()
    {
        // The spectra come from the A549 database, so searching it as a transient database
        // gives transient identifications and transient protein groups.
        string outputFolder = FreshOutputFolder("ParallelSearchAllOutputs");
        try
        {
            string a549 = Path.GetFileNameWithoutExtension(BaseDatabasePath);
            var task = new TaskLayer.ParallelSearch.ParallelSearchTask(
            [
                new DbForTask(BaseDatabasePath, false),
                new DbForTask(TransientDatabasePaths[1], false)
            ]);
            task.ParallelSearchParameters.WriteTransientResultsOnly = false;
            task.ParallelSearchParameters.WriteSpectralLibrary = true;
            task.ParallelSearchParameters.WriteTransientSpectralLibrary = true;
            task.ParallelSearchParameters.CompressTransientSearchOutputs = true;

            task.RunTask(outputFolder, [new DbForTask(TransientDatabasePaths[0], false)], [SpectraPath], "parallel");

            var summary = ReadSummary(outputFolder);
            string databaseFolder = Path.Combine(outputFolder, a549);
            string[] writtenFiles = Directory.GetFiles(databaseFolder).Select(Path.GetFileName).ToArray();

            Assert.Multiple(() =>
            {
                Assert.That(summary.Keys, Has.Count.EqualTo(2));
                Assert.That(summary[a549].TargetPsmsFromTransientDbAtQValueThreshold, Is.GreaterThan(0));
                Assert.That(summary[a549].TargetProteinGroupsFromTransientDb, Is.GreaterThan(0));
                Assert.That(task.TotalMs2Scans, Is.GreaterThan(0));

                // every output is compressed in place
                Assert.That(writtenFiles, Has.Member($"{a549}_AllPSMs.psmtsv.gz"));
                Assert.That(writtenFiles, Has.Member("AllPSMs.psmtsv.gz"));
                Assert.That(writtenFiles, Has.Member($"{a549}_AllPeptides.psmtsv.gz"));
                Assert.That(writtenFiles, Has.Member("AllPeptides.psmtsv.gz"));
                Assert.That(writtenFiles, Has.Member($"{a549}_AllProteinGroups.tsv.gz"));
                Assert.That(writtenFiles, Has.Member("AllProteinGroups.tsv.gz"));
                Assert.That(writtenFiles, Has.Member($"{a549}_SpectralLibrary.msp.gz"));
                Assert.That(writtenFiles, Has.Member($"AllPeptidesAnd_{a549}_SpectralLibrary.msp.gz"));
                Assert.That(writtenFiles, Has.Member("results.txt.gz"));
            });
        }
        finally
        {
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);
        }
    }

    [Test]
    public void RunTask_ResumeSearchesOnlyNewDatabases_AndOverwriteRedoesExistingOnes()
    {
        string outputFolder = FreshOutputFolder("ParallelSearchResume");
        try
        {
            var baseDatabases = new List<DbForTask> { new(BaseDatabasePath, false) };
            var spectra = new List<string> { SpectraPath };
            string first = Path.GetFileNameWithoutExtension(TransientDatabasePaths[0]);
            string second = Path.GetFileNameWithoutExtension(TransientDatabasePaths[1]);
            string firstPsmFile = Path.Combine(outputFolder, first, $"{first}_AllPSMs.psmtsv");

            // 1: only the first database
            var onlyFirst = new TaskLayer.ParallelSearch.ParallelSearchTask([new DbForTask(TransientDatabasePaths[0], false)]);
            onlyFirst.RunTask(outputFolder, baseDatabases, spectra, "parallel");
            var firstWriteTime = File.GetLastWriteTimeUtc(firstPsmFile);
            Assert.That(ReadSummary(outputFolder).Keys, Is.EquivalentTo(new[] { first }));

            // 2: both databases, keeping what exists: the first is skipped and the second is searched
            CreateTask(overwrite: false).RunTask(outputFolder, baseDatabases, spectra, "parallel");
            Assert.Multiple(() =>
            {
                Assert.That(ReadSummary(outputFolder).Keys, Is.EquivalentTo(new[] { first, second }));
                Assert.That(File.GetLastWriteTimeUtc(firstPsmFile), Is.EqualTo(firstWriteTime));
                Assert.That(File.Exists(Path.Combine(outputFolder, second, "results.txt")), Is.True);
            });

            // 3: the first database again, overwriting: its folder is rebuilt and the second row is kept
            string marker = Path.Combine(outputFolder, first, "left-over.txt");
            File.WriteAllText(marker, "stale");
            var overwriteFirst = new TaskLayer.ParallelSearch.ParallelSearchTask([new DbForTask(TransientDatabasePaths[0], false)]);
            overwriteFirst.ParallelSearchParameters.OverwriteTransientSearchOutputs = true;
            overwriteFirst.RunTask(outputFolder, baseDatabases, spectra, "parallel");

            string[] summaryLines = File.ReadAllLines(
                Path.Combine(outputFolder, TransientDatabaseResultsManager.SummaryResultsFileName));
            Assert.Multiple(() =>
            {
                Assert.That(File.Exists(marker), Is.False);
                Assert.That(File.Exists(firstPsmFile), Is.True);
                Assert.That(summaryLines.Length, Is.EqualTo(3));
                Assert.That(ReadSummary(outputFolder).Keys, Is.EquivalentTo(new[] { first, second }));
            });
        }
        finally
        {
            if (Directory.Exists(outputFolder))
                Directory.Delete(outputFolder, true);
        }
    }

    /// <summary>
    /// The PSM and peptide files follow a normal search's write settings: with decoys and
    /// high q-value matches turned off, neither appears in the baseline or per-database files.
    /// </summary>
    [Test]
    public void RunTask_PsmAndPeptideFilesHonourTheWriteSettings()
    {
        string outputFolder = FreshOutputFolder("ParallelSearchWriteSettings");
        try
        {
            string a549 = Path.GetFileNameWithoutExtension(BaseDatabasePath);
            var task = new TaskLayer.ParallelSearch.ParallelSearchTask([new DbForTask(BaseDatabasePath, false)]);
            task.ParallelSearchParameters.WriteTransientResultsOnly = false;
            task.ParallelSearchParameters.WriteDecoys = false;
            task.ParallelSearchParameters.WriteHighQValuePsms = false;

            task.RunTask(outputFolder, [new DbForTask(TransientDatabasePaths[0], false)], [SpectraPath], "parallel");

            string[] files =
            [
                Path.Combine(outputFolder, "BasePSMs.psmtsv"),
                Path.Combine(outputFolder, "BasePeptides.psmtsv"),
                Path.Combine(outputFolder, a549, $"{a549}_AllPSMs.psmtsv"),
                Path.Combine(outputFolder, a549, "AllPSMs.psmtsv"),
                Path.Combine(outputFolder, a549, $"{a549}_AllPeptides.psmtsv"),
                Path.Combine(outputFolder, a549, "AllPeptides.psmtsv"),
            ];

            Assert.Multiple(() =>
            {
                foreach (string file in files)
                {
                    string[] lines = File.ReadAllLines(file);
                    string[] header = lines[0].Split('\t');
                    int decoyColumn = System.Array.IndexOf(header, "Decoy/Contaminant/Target");
                    int qValueColumn = System.Array.IndexOf(header, "QValue");
                    var rows = lines.Skip(1).Select(l => l.Split('\t')).ToList();

                    Assert.That(rows, Is.Not.Empty, file);
                    Assert.That(rows.Where(r => r[decoyColumn].Contains('D')), Is.Empty, $"decoys written to {file}");
                    Assert.That(rows.Select(r => double.Parse(r[qValueColumn], System.Globalization.CultureInfo.InvariantCulture)),
                        Has.All.LessThanOrEqualTo(0.01), $"high q-value matches written to {file}");
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
    public void RunTask_ReportsRetentionTimeErrorsInMinutes()
    {
        // The A549 spectra searched against A549 as both base and transient database: the base search gives the
        // anchors for the per-file Chronologer calibration, the transient search gives identifications to score.
        string outputFolder = FreshOutputFolder("ParallelSearchRetentionTime");
        try
        {
            string a549 = Path.GetFileNameWithoutExtension(BaseDatabasePath);
            var task = new TaskLayer.ParallelSearch.ParallelSearchTask([new DbForTask(BaseDatabasePath, false)]);

            task.RunTask(outputFolder, [new DbForTask(BaseDatabasePath, false)], [SpectraPath], "parallel");

            string summaryPath = Path.Combine(outputFolder, TransientDatabaseResultsManager.SummaryResultsFileName);
            string[] header = File.ReadAllLines(summaryPath)[0].Split(',');
            var row = ReadSummary(outputFolder)[a549];
            TestContext.Out.WriteLine($"PSM errors: n={row.Psm_AllRtErrors.Length} mean |error|={row.Psm_MeanAbsoluteRtError} " +
                $"min={row.Psm_AllRtErrors.DefaultIfEmpty(double.NaN).Min()} max={row.Psm_AllRtErrors.DefaultIfEmpty(double.NaN).Max()}");
            TestContext.Out.WriteLine($"peptide errors: n={row.Peptide_AllRtErrors.Length} mean |error|={row.Peptide_MeanAbsoluteRtError}");

            Assert.Multiple(() =>
            {
                Assert.That(header, Does.Contain("Psm_MeanAbsoluteRtError"));
                Assert.That(header, Does.Not.Contain("Psm_MeanAbsoluteRtError"));
                Assert.That(row.Psm_AllRtErrors, Is.Not.Empty);
                Assert.That(row.Peptide_AllRtErrors, Is.Not.Empty);
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
