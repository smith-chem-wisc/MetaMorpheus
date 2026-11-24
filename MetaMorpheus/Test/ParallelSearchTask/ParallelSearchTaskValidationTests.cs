using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using NUnit.Framework;
using TaskLayer;
using TaskLayer.ParallelSearch;
using Transcriptomics.Digestion;

namespace Test.ParallelSearchTask;

/// <summary>
/// Settings the parallel search cannot honour must be rejected before anything is loaded or written.
/// </summary>
[TestFixture]
public class ParallelSearchTaskValidationTests
{
    private static readonly string TestDataFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
    private static readonly string SpectraPath = Path.Combine(TestDataFolder, "TaGe_SA_A549_3_snip.mzML");
    private static readonly string BaseDatabasePath = Path.Combine(TestDataFolder, "TaGe_SA_A549_3_snip.fasta");
    private static readonly string TransientDatabasePath = Path.Combine(TestDataFolder, "smalldb.fasta");
    private static readonly string ValidationFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "ParallelSearchValidation");

    private static TaskLayer.ParallelSearch.ParallelSearchTask CreateTask(CommonParameters? commonParameters = null)
    {
        var task = new TaskLayer.ParallelSearch.ParallelSearchTask([new DbForTask(TransientDatabasePath, false)]);
        if (commonParameters != null)
        {
            task.CommonParameters = commonParameters;
        }
        return task;
    }

    /// <summary>
    /// Runs the task, expects a <see cref="MetaMorpheusException"/>, and checks that the task wrote nothing of its
    /// own: RunTask always leaves results.txt (with the error) in the output folder, so that is the only file allowed.
    /// </summary>
    private static MetaMorpheusException RunAndExpectRejection(TaskLayer.ParallelSearch.ParallelSearchTask task, string caseName,
        List<string>? spectraPaths = null)
    {
        string outputFolder = Path.Combine(ValidationFolder, caseName, "Output");
        if (Directory.Exists(outputFolder))
        {
            Directory.Delete(outputFolder, true);
        }
        Directory.CreateDirectory(outputFolder);

        var exception = Assert.Throws<MetaMorpheusException>(() => task.RunTask(outputFolder,
            [new DbForTask(BaseDatabasePath, false)], spectraPaths ?? [SpectraPath], "parallel"))!;

        string[] written = Directory.GetFileSystemEntries(outputFolder).Select(Path.GetFileName).ToArray()!;
        Assert.That(written, Is.EquivalentTo(new[] { "results.txt" }), "the task wrote output before rejecting its settings");
        return exception;
    }

    [OneTimeTearDown]
    public void CleanUp()
    {
        if (Directory.Exists(ValidationFolder))
        {
            Directory.Delete(ValidationFolder, true);
        }
    }

    [Test]
    public void RunTask_WithNoTransientDatabases_IsRejected()
    {
        var task = new TaskLayer.ParallelSearch.ParallelSearchTask();

        var exception = RunAndExpectRejection(task, nameof(RunTask_WithNoTransientDatabases_IsRejected));

        Assert.That(exception.Message, Does.Contain("at least one transient database"));
    }

    [Test]
    public void RunTask_WithTransientDatabasesSharingAFileName_IsRejectedAndNamesThem()
    {
        string otherFolder = Path.Combine(ValidationFolder, "OtherFolder");
        Directory.CreateDirectory(otherFolder);
        string sameNamePath = Path.Combine(otherFolder, "smalldb.xml");
        var task = CreateTask();
        task.ParallelSearchParameters.TransientDatabases.Add(new DbForTask(sameNamePath, false));

        var exception = RunAndExpectRejection(task, nameof(RunTask_WithTransientDatabasesSharingAFileName_IsRejectedAndNamesThem));

        Assert.Multiple(() =>
        {
            Assert.That(exception.Message, Does.Contain("different file name"));
            Assert.That(exception.Message, Does.Contain("smalldb"));
            Assert.That(exception.Message, Does.Contain(sameNamePath));
        });
    }

    [Test]
    public void RunTask_InRnaMode_IsRejected()
    {
        try
        {
            var task = CreateTask(new CommonParameters(digestionParams: new RnaDigestionParams("RNase T1")));

            var exception = RunAndExpectRejection(task, nameof(RunTask_InRnaMode_IsRejected));

            Assert.That(exception.Message, Does.Contain("RNA"));
        }
        finally
        {
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
        }
    }

    [Test]
    [TestCase(DissociationType.Autodetect, "autodetect")]
    [TestCase(DissociationType.LowCID, "Low CID")]
    public void RunTask_WithUnsupportedDissociationType_IsRejected(DissociationType dissociationType, string expectedText)
    {
        var task = CreateTask(new CommonParameters(dissociationType: dissociationType));

        var exception = RunAndExpectRejection(task, $"{nameof(RunTask_WithUnsupportedDissociationType_IsRejected)}_{dissociationType}");

        Assert.That(exception.Message, Does.Contain(expectedText));
    }

    [Test]
    [TestCase(DissociationType.Autodetect, "autodetect")]
    [TestCase(DissociationType.LowCID, "Low CID")]
    public void RunTask_WithUnsupportedFileSpecificDissociationType_IsRejected(DissociationType dissociationType, string expectedText)
    {
        // A file-specific toml sits beside the spectra file it applies to, so work on a copy
        string spectraFolder = Path.Combine(ValidationFolder, $"FileSpecific_{dissociationType}");
        Directory.CreateDirectory(spectraFolder);
        string spectraCopy = Path.Combine(spectraFolder, Path.GetFileName(SpectraPath));
        File.Copy(SpectraPath, spectraCopy, true);
        File.WriteAllText(Path.ChangeExtension(spectraCopy, ".toml"), $"DissociationType = \"{dissociationType}\"{Environment.NewLine}");

        var exception = RunAndExpectRejection(CreateTask(),
            $"{nameof(RunTask_WithUnsupportedFileSpecificDissociationType_IsRejected)}_{dissociationType}", [spectraCopy]);

        Assert.Multiple(() =>
        {
            Assert.That(exception.Message, Does.Contain(expectedText));
            Assert.That(exception.Message, Does.Contain("file-specific"));
        });
    }

    [Test]
    public void RunTask_WithFileSpecificDigestionOrFragmentationThatDiffersFromTheTask_IsRejectedAndNamesThem()
    {
        // the transient engine digests and fragments with the task's settings, so a file asking for ETD or
        // different missed cleavages would be searched as HCD with the task's digestion and no warning
        string spectraFolder = Path.Combine(ValidationFolder, "FileSpecificDigestion");
        Directory.CreateDirectory(spectraFolder);
        string spectraCopy = Path.Combine(spectraFolder, Path.GetFileName(SpectraPath));
        File.Copy(SpectraPath, spectraCopy, true);
        File.WriteAllText(Path.ChangeExtension(spectraCopy, ".toml"),
            $"DissociationType = \"ETD\"{Environment.NewLine}MaxMissedCleavages = 5{Environment.NewLine}");

        var task = CreateTask(new CommonParameters(dissociationType: DissociationType.HCD));
        var exception = RunAndExpectRejection(task,
            nameof(RunTask_WithFileSpecificDigestionOrFragmentationThatDiffersFromTheTask_IsRejectedAndNamesThem), [spectraCopy]);

        Assert.Multiple(() =>
        {
            Assert.That(exception.Message, Does.Contain(Path.GetFileName(spectraCopy)));
            Assert.That(exception.Message, Does.Contain("DissociationType (ETD for the file, HCD for the task)"));
            Assert.That(exception.Message, Does.Contain("MaxMissedCleavages (5 for the file, 2 for the task)"));
        });
    }

    [Test]
    public void RunTask_WithMostAbundantMassDiffAcceptor_IsRejected()
    {
        var task = CreateTask();
        task.ParallelSearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.MostAbundant_Exact;

        var exception = RunAndExpectRejection(task, nameof(RunTask_WithMostAbundantMassDiffAcceptor_IsRejected));

        Assert.That(exception.Message, Does.Contain("most-abundant"));
    }

    [Test]
    [TestCase(SearchType.Modern)]
    [TestCase(SearchType.NonSpecific)]
    public void RunTask_WithNonClassicSearchType_IsRejected(SearchType searchType)
    {
        var task = CreateTask();
        task.ParallelSearchParameters.SearchType = searchType;

        var exception = RunAndExpectRejection(task, $"{nameof(RunTask_WithNonClassicSearchType_IsRejected)}_{searchType}");

        Assert.That(exception.Message, Does.Contain("classic search").And.Contain(searchType.ToString()));
    }

    [Test]
    public void RunTask_WithMostAbundantPrecursorMatchMode_IsRejected()
    {
        var task = CreateTask(new CommonParameters(precursorMassMatchMode: PrecursorMassMatchMode.MostAbundant));

        var exception = RunAndExpectRejection(task, nameof(RunTask_WithMostAbundantPrecursorMatchMode_IsRejected));

        Assert.That(exception.Message, Does.Contain("most-abundant"));
    }

    [Test]
    public void SearchParametersSetter_WithPlainSearchParameters_KeepsTheAssignedSettings()
    {
        var task = new TaskLayer.ParallelSearch.ParallelSearchTask();

        task.SearchParameters = new SearchParameters
        {
            SearchTarget = false,
            DecoyType = UsefulProteomicsDatabases.DecoyType.None,
            WriteDecoys = false,
            ModPeptidesAreDifferent = true,
        };

        Assert.Multiple(() =>
        {
            Assert.That(task.SearchParameters, Is.SameAs(task.ParallelSearchParameters));
            Assert.That(task.ParallelSearchParameters.SearchTarget, Is.False);
            Assert.That(task.ParallelSearchParameters.DecoyType, Is.EqualTo(UsefulProteomicsDatabases.DecoyType.None));
            Assert.That(task.ParallelSearchParameters.WriteDecoys, Is.False);
            Assert.That(task.ParallelSearchParameters.ModPeptidesAreDifferent, Is.True);
        });
    }
}
