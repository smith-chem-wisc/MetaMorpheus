using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Deconvolution;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test;

[TestFixture]
public static class FromFileDeconvolutionSearchTests
{
    [Test]
    public static void SearchTask_UsesFlashDeconvFeaturesForThreeRawFiles()
    {
        string testDirectory = TestContext.CurrentContext.TestDirectory;
        string dataDirectory = Path.Combine(testDirectory, "TestData");
        string featureDirectory = Path.Combine(dataDirectory, "file-specific-decon");
        string[] inputBaseNames =
        [
            "TaGe_SA_A549_3_snip",
            "TaGe_SA_A549_3_snip_2",
            "TaGe_SA_HeLa_04_subset_longestSeq"
        ];

        var rawFiles = inputBaseNames
            .Select(name => Path.Combine(dataDirectory, name + ".mzML"))
            .ToList();
        var map = new SearchFeatureFileMap(inputBaseNames.Select((name, index) =>
            new SearchFeatureFileMapEntry(
                massSpecFilePath: rawFiles[index],
                featureFilePath: Path.Combine(featureDirectory, name + "_ms1.feature"))));
        var mappedParameters = new FeatureMappedFromFileDeconvolutionParameters(map, minCharge: 1, maxCharge: 20)
        {
            UseGenericScore = true
        };

        for (int i = 0; i < rawFiles.Count; i++)
        {
            string rawFilePath = rawFiles[i];
            Assert.That(map.TryGetFeaturePathForMassSpecFile(rawFilePath, out var featureFilePath), Is.True);
            Assert.That(File.Exists(rawFilePath), Is.True, rawFilePath);
            Assert.That(File.Exists(featureFilePath), Is.True, featureFilePath);
            var fileParameters = (FromFileDeconvolutionParameters)
                mappedParameters.ToDeconvolutionParameters(rawFilePath);
            Assert.That(fileParameters.Features, Is.Not.Empty, featureFilePath);
        }

        var searchTask = new SearchTask
        {
            SearchParameters = new SearchParameters { Normalize = false },
            CommonParameters = new CommonParameters(
                useProvidedPrecursorInfo: false,
                precursorDeconParams: mappedParameters)
        };
        string databasePath = Path.Combine(dataDirectory, "TaGe_SA_A549_3_snip.fasta");
        string temporaryRoot = Path.Combine(testDirectory, "FromFileDeconvolutionSearch_");
        string firstRunRoot = Path.Combine(temporaryRoot, "Manual");
        string secondRunRoot = Path.Combine(temporaryRoot, "TomlRoundTrip");
        string firstOutputDirectory = CreateRunDirectories(firstRunRoot);
        string secondOutputDirectory = CreateRunDirectories(secondRunRoot);

        try
        {
            searchTask.RunTask(
                firstOutputDirectory,
                new List<DbForTask> { new(databasePath, false) },
                rawFiles,
                "FromFileDeconvolutionSearch");

            string taskTomlPath = Path.Combine(
                firstRunRoot,
                "Task Settings",
                "FromFileDeconvolutionSearchconfig.toml");
            Assert.That(File.Exists(taskTomlPath), Is.True, taskTomlPath);
            string taskToml = File.ReadAllText(taskTomlPath);
            Assert.That(taskToml, Does.Contain("FeatureFileMap]"));
            Assert.That(taskToml, Does.Contain("'TaGe_SA_A549_3_snip.mzML' = "));
            Assert.That(taskToml, Does.Contain("'TaGe_SA_A549_3_snip_2.mzML' = "));
            Assert.That(taskToml, Does.Contain("'TaGe_SA_HeLa_04_subset_longestSeq.mzML' = "));
            Assert.That(taskToml, Does.Not.Contain("Entries"));
            Assert.That(taskToml, Does.Contain("file-specific-decon"));
            Assert.That(taskToml, Does.Not.Contain(testDirectory));

            var roundTrippedTask = Toml.ReadString<SearchTask>(taskToml, MetaMorpheusTask.tomlConfig);
            var roundTrippedParameters = roundTrippedTask.CommonParameters.PrecursorDeconvolutionParameters
                as FeatureMappedFromFileDeconvolutionParameters;
            Assert.That(roundTrippedParameters, Is.Not.Null);
            Assert.That(roundTrippedParameters.FeatureFileMap, Is.EqualTo(map));
            Assert.That(roundTrippedParameters.UseGenericScore, Is.EqualTo(mappedParameters.UseGenericScore));

            roundTrippedTask.RunTask(
                secondOutputDirectory,
                new List<DbForTask> { new(databasePath, false) },
                rawFiles,
                "FromFileDeconvolutionSearch");

            string firstPsmsPath = Path.Combine(firstOutputDirectory, "AllPSMs.psmtsv");
            string secondPsmsPath = Path.Combine(secondOutputDirectory, "AllPSMs.psmtsv");
            Assert.That(File.Exists(firstPsmsPath), Is.True);
            Assert.That(File.Exists(secondPsmsPath), Is.True);
            byte[] firstPsms = File.ReadAllBytes(firstPsmsPath);
            byte[] secondPsms = File.ReadAllBytes(secondPsmsPath);
            Assert.That(firstPsms.Length, Is.GreaterThan(0));
            Assert.That(firstPsms.SequenceEqual(secondPsms), Is.True, "TOML round-trip search changed the PSM file bytes.");
        }
        catch (Exception ex)
        {
            Assert.Fail($"Search task failed with exception: {ex}");
        }
        finally
        {
            if (Directory.Exists(temporaryRoot))
                Directory.Delete(temporaryRoot, true);
        }
    }

    /// <summary>
    /// End-to-end coverage for the file-specific route: each raw file carries its own hand-written
    /// companion TOML beside it (never next to the repository fixtures), containing native mzLib
    /// <see cref="FromFileDeconvolutionParameters"/>. The search must consume those per-raw feature
    /// files and produce PSM output without falling back to the task-wide precursor deconvolution.
    /// </summary>
    [Test, NonParallelizable] // subscribes to the process-wide MetaMorpheusTask.WarnHandler
    public static void SearchTask_FileSpecificFromFile_UsesFlashDeconvFeaturesForThreeRawFiles()
    {
        string temporaryRoot = CreateIsolatedTemporaryRoot("FileSpecificSuccess");
        string inputDirectory = Path.Combine(temporaryRoot, "Input");
        string outputDirectory = CreateRunDirectories(Path.Combine(temporaryRoot, "FileSpecificRun"));

        try
        {
            List<string> rawFiles = CopyApprovedRawFixtures(inputDirectory);

            // Standalone file-specific route: the task-wide precursor deconvolution is Classic, so a
            // per-raw fallback would show up in the captured warnings and in the resolved parameters.
            var commonParameters = CreateCommonParameters(
                new ClassicDeconvolutionParameters(MinCharge, MaxCharge, 4, 3));
            for (int i = 0; i < rawFiles.Count; i++)
            {
                string featurePath = FeatureFixturePath(InputBaseNames[i]);
                Assert.That(File.Exists(featurePath), Is.True, $"Approved feature fixture missing: {featurePath}");
                WriteCompanionToml(rawFiles[i], featurePath);
                AssertCompanionOverridesTaskWide(
                    LoadCompanionSettings(rawFiles[i]), commonParameters, rawFiles[i], featurePath);
            }

            List<string> warnings = RunSearchCapturingWarnings(
                CreateSearchTask(commonParameters), outputDirectory, rawFiles);

            AssertNoFileSpecificFallbackWarnings(warnings);

            string psmsPath = Path.Combine(outputDirectory, "AllPSMs.psmtsv");
            Assert.That(File.Exists(psmsPath), Is.True, psmsPath);
            Assert.That(CountPsmRows(psmsPath), Is.GreaterThan(0),
                "The file-specific FromFile search must write at least one PSM row.");
        }
        finally
        {
            if (Directory.Exists(temporaryRoot))
                Directory.Delete(temporaryRoot, true);
        }
    }

    /// <summary>
    /// Byte-equivalence coverage: the same raw files, FASTA, raw order, CommonParameters (including
    /// the task-wide embedded raw-to-feature map), SearchParameters, feature pairings, charge bounds,
    /// polarity and <c>UseGenericScore</c> are run once through the existing task-map route and once
    /// through per-raw companion TOMLs with the equivalent native FromFile parameters. Only the two
    /// <c>AllPSMs.psmtsv</c> files are compared, byte for byte.
    /// </summary>
    [Test, NonParallelizable] // subscribes to the process-wide MetaMorpheusTask.WarnHandler
    public static void SearchTask_FileSpecificFromFile_MatchesTaskMapOutputBytes()
    {
        string temporaryRoot = CreateIsolatedTemporaryRoot("EquivalentRoutes");
        string inputDirectory = Path.Combine(temporaryRoot, "Input");
        string mappedOutputDirectory = CreateRunDirectories(Path.Combine(temporaryRoot, "Mapped"));
        string fileSpecificOutputDirectory = CreateRunDirectories(Path.Combine(temporaryRoot, "FileSpecific"));

        try
        {
            List<string> rawFiles = CopyApprovedRawFixtures(inputDirectory);
            var map = CreateFixtureMap(rawFiles);
            var mappedParameters = new FeatureMappedFromFileDeconvolutionParameters(
                map, minCharge: MinCharge, maxCharge: MaxCharge, polarity: Polarity.Positive)
            {
                UseGenericScore = true
            };

            var commonParameters = CreateCommonParameters(mappedParameters);

            List<string> mappedWarnings = RunSearchCapturingWarnings(
                CreateSearchTask(commonParameters), mappedOutputDirectory, rawFiles);
            AssertNoFileSpecificFallbackWarnings(mappedWarnings);

            // Write companion TOMLs only after the mapped run so the mapped route cannot consume them.
            for (int i = 0; i < rawFiles.Count; i++)
            {
                Assert.That(File.Exists(Path.Combine(inputDirectory, InputBaseNames[i] + ".toml")), Is.False,
                    "Companion TOMLs must not exist while the mapped route runs");
                WriteCompanionToml(rawFiles[i], FeatureFixturePath(InputBaseNames[i]));
                AssertCompanionOverridesTaskWide(
                    LoadCompanionSettings(rawFiles[i]), commonParameters, rawFiles[i], FeatureFixturePath(InputBaseNames[i]));
            }

            List<string> fileSpecificWarnings = RunSearchCapturingWarnings(
                CreateSearchTask(commonParameters), fileSpecificOutputDirectory, rawFiles);
            AssertNoFileSpecificFallbackWarnings(fileSpecificWarnings);

            string mappedPsmsPath = Path.Combine(mappedOutputDirectory, "AllPSMs.psmtsv");
            string fileSpecificPsmsPath = Path.Combine(fileSpecificOutputDirectory, "AllPSMs.psmtsv");
            Assert.That(File.Exists(mappedPsmsPath), Is.True, mappedPsmsPath);
            Assert.That(File.Exists(fileSpecificPsmsPath), Is.True, fileSpecificPsmsPath);

            byte[] mappedPsms = File.ReadAllBytes(mappedPsmsPath);
            byte[] fileSpecificPsms = File.ReadAllBytes(fileSpecificPsmsPath);
            Assert.That(mappedPsms.Length, Is.GreaterThan(0));
            Assert.That(fileSpecificPsms, Is.EqualTo(mappedPsms),
                "The per-file companion TOML route changed the AllPSMs.psmtsv bytes relative to the task-map route.");
        }
        finally
        {
            if (Directory.Exists(temporaryRoot))
                Directory.Delete(temporaryRoot, true);
        }
    }

    private const int MinCharge = 1;
    private const int MaxCharge = 20;
    private const string DisplayName = "FromFileDeconvolutionSearch";
    private const string FileSpecificFallbackWarningPrefix = "File-specific precursor deconvolution";

    /// <summary>The three approved raw/FlashDeconv fixture pairs, kept in a fixed order.</summary>
    private static readonly string[] InputBaseNames =
    [
        "TaGe_SA_A549_3_snip",
        "TaGe_SA_A549_3_snip_2",
        "TaGe_SA_HeLa_04_subset_longestSeq"
    ];

    private static string TestDataDirectory =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");

    private static string FeatureFixturePath(string inputBaseName) =>
        Path.Combine(TestDataDirectory, "file-specific-decon", inputBaseName + "_ms1.feature");

    /// <summary>
    /// Isolated root under NUnit's test directory; the raw fixtures are copied into it because the
    /// companion <c>&lt;raw basename&gt;.toml</c> files must sit beside the raw files and must never be
    /// written next to the repository fixtures.
    /// </summary>
    private static string CreateIsolatedTemporaryRoot(string scenario)
    {
        string root = Path.Combine(
            TestContext.CurrentContext.TestDirectory,
            $"FromFileDeconvolutionSearch_{scenario}_{Guid.NewGuid():N}");
        Directory.CreateDirectory(root);
        return root;
    }

    private static List<string> CopyApprovedRawFixtures(string inputDirectory)
    {
        Directory.CreateDirectory(inputDirectory);
        var copiedRawFiles = new List<string>();
        foreach (string inputBaseName in InputBaseNames)
        {
            string source = Path.Combine(TestDataDirectory, inputBaseName + ".mzML");
            Assert.That(File.Exists(source), Is.True, $"Approved raw fixture missing: {source}");
            string destination = Path.Combine(inputDirectory, inputBaseName + ".mzML");
            File.Copy(source, destination);
            copiedRawFiles.Add(destination);
        }

        return copiedRawFiles;
    }

    private static List<DbForTask> CreateDatabaseList() =>
        new() { new(Path.Combine(TestDataDirectory, "TaGe_SA_A549_3_snip.fasta"), false) };

    private static SearchFeatureFileMap CreateFixtureMap(IReadOnlyList<string> copiedRawFiles) =>
        new()
        {
            Entries = copiedRawFiles
                .Select((rawFile, index) => new SearchFeatureFileMapEntry(rawFile, FeatureFixturePath(InputBaseNames[index])))
                .ToList()
        };

    private static SearchTask CreateSearchTask(CommonParameters commonParameters) =>
        new()
        {
            SearchParameters = new SearchParameters { Normalize = false },
            CommonParameters = commonParameters
        };

    private static CommonParameters CreateCommonParameters(DeconvolutionParameters taskWidePrecursorDeconvolution) =>
        new(
            doPrecursorDeconvolution: true,
            useProvidedPrecursorInfo: false,
            precursorDeconParams: taskWidePrecursorDeconvolution);

    /// <summary>
    /// Writes the production-shaped companion TOML next to the (copied) raw file: a
    /// <see cref="FileSpecificParameters"/> whose only value is a native mzLib FromFile precursor
    /// deconvolution pointing at the absolute feature-fixture path.
    /// </summary>
    private static void WriteCompanionToml(
        string rawFilePath, string featureFilePath, int minCharge = MinCharge, int maxCharge = MaxCharge)
    {
        var fileSpecificParameters = new FileSpecificParameters
        {
            PrecursorDeconvolutionParameters = new FromFileDeconvolutionParameters(
                featureFilePath, minCharge, maxCharge, Polarity.Positive)
            {
                UseGenericScore = true
            }
        };
        File.WriteAllText(
            Path.ChangeExtension(rawFilePath, ".toml"),
            Toml.WriteString(fileSpecificParameters, MetaMorpheusTask.tomlConfig));
    }

    private static FileSpecificParameters LoadCompanionSettings(string rawFilePath) =>
        new(Toml.ReadFile(Path.ChangeExtension(rawFilePath, ".toml"), MetaMorpheusTask.tomlConfig));

    /// <summary>
    /// Resolves the companion settings through the same canonical resolver the search uses and
    /// asserts the native file-specific FromFile object wins, carries the expected values and has
    /// preloaded a non-empty feature collection.
    /// </summary>
    private static void AssertCompanionOverridesTaskWide(
        FileSpecificParameters companionSettings, CommonParameters commonParameters,
        string rawFilePath, string expectedFeaturePath)
    {
        Assert.That(companionSettings.PrecursorDeconvolutionParameters, Is.TypeOf<FromFileDeconvolutionParameters>());
        var native = (FromFileDeconvolutionParameters)companionSettings.PrecursorDeconvolutionParameters;
        Assert.That(native.FilePath, Is.EqualTo(expectedFeaturePath));

        var resolved = MetaMorpheusTask.SetAllFileSpecificCommonParams(commonParameters, companionSettings, rawFilePath);
        Assert.That(resolved.PrecursorDeconvolutionParameters, Is.SameAs(companionSettings.PrecursorDeconvolutionParameters),
            $"The companion TOML for '{rawFilePath}' must win over the task-wide parameters");
        var resolvedNative = (FromFileDeconvolutionParameters)resolved.PrecursorDeconvolutionParameters;
        Assert.That(resolvedNative.MinAssumedChargeState, Is.EqualTo(MinCharge));
        Assert.That(resolvedNative.MaxAssumedChargeState, Is.EqualTo(MaxCharge));
        Assert.That(resolvedNative.Polarity, Is.EqualTo(Polarity.Positive));
        Assert.That(resolvedNative.UseGenericScore, Is.True);
        Assert.That(resolvedNative.Features, Is.Not.Empty,
            $"The preflight must load MS1 features for '{expectedFeaturePath}'");
    }

    private static List<string> RunSearchCapturingWarnings(
        SearchTask searchTask, string outputDirectory, List<string> rawFiles)
    {
        var warnings = new List<string>();
        EventHandler<StringEventArgs> handler = (_, e) => warnings.Add(e.S);
        MetaMorpheusTask.WarnHandler += handler;
        try
        {
            searchTask.RunTask(outputDirectory, CreateDatabaseList(), rawFiles, DisplayName);
        }
        finally
        {
            MetaMorpheusTask.WarnHandler -= handler;
        }

        return warnings;
    }

    private static void AssertNoFileSpecificFallbackWarnings(List<string> warnings) =>
        Assert.That(
            warnings.Where(warning => warning.Contains(FileSpecificFallbackWarningPrefix)),
            Is.Empty,
            "A valid companion TOML must not fall back to the task-wide precursor deconvolution");

    private static int CountPsmRows(string psmTsvPath) =>
        File.ReadLines(psmTsvPath).Skip(1).Count(line => !string.IsNullOrWhiteSpace(line));

    private static string CreateRunDirectories(string runRoot)
    {
        string outputDirectory = Path.Combine(runRoot, "Output");
        Directory.CreateDirectory(outputDirectory);
        Directory.CreateDirectory(Path.Combine(runRoot, "Task Settings"));
        return outputDirectory;
    }
}
