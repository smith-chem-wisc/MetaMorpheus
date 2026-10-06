using EngineLayer;
using EngineLayer.Deconvolution;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test;

/// <summary>
/// Runtime preflight/fallback coverage for a file-specific native mzLib
/// <see cref="FromFileDeconvolutionParameters"/> precursor override.
/// <para>
/// <see cref="MetaMorpheusTask.SetAllFileSpecificCommonParams"/> must force the lazy mzLib feature
/// load before any MS2 scan processing, keep a valid override for its raw file only, and isolate a
/// bad override to its own raw with a single warning naming the raw path, the selected feature path,
/// the failure reason and the task-wide fallback that was used instead. The companion TOML must
/// never be rewritten by the runtime fallback.
/// </para>
/// </summary>
[TestFixture, NonParallelizable] // subscribes to the static MetaMorpheusTask.WarnHandler
public static class FileSpecificFromFileFallbackTests
{
    private const string RawWithOverride = @"E:\data\raw_with_override.mzML";
    private const string RawWithoutOverride = @"E:\data\raw_without_override.mzML";

    private static string TestDataDirectory =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "file-specific-decon");

    private static string ValidMs1FeaturePath =>
        Path.Combine(TestDataDirectory, "TaGe_SA_A549_3_snip_ms1.feature");

    private static string SecondValidMs1FeaturePath =>
        Path.Combine(TestDataDirectory, "TaGe_SA_A549_3_snip_2_ms1.feature");

    private static string NonMs1FeaturePath =>
        Path.Combine(TestDataDirectory, "TaGe_SA_A549_3_snip_ms2.feature");

    private static string NonexistentFeaturePath =>
        Path.Combine(TestDataDirectory, "does_not_exist_ms1.feature");

    [Test]
    public static void FileSpecificFromFile_Fallback_MissingFeatureFile_WarnsAndFallsBack()
    {
        var commonParameters = ClassicCommonParameters();
        var fileSpecific = FromFileSettings(NonexistentFeaturePath);

        List<string> warnings = CaptureWarnings(() =>
        {
            var resolved = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                commonParameters, fileSpecific, RawWithOverride);

            Assert.That(resolved.PrecursorDeconvolutionParameters, Is.TypeOf<ClassicDeconvolutionParameters>());
            Assert.That(resolved.PrecursorDeconvolutionParameters, Is.SameAs(commonParameters.PrecursorDeconvolutionParameters),
                "A failed file-specific override must fall back to the task-wide precursor parameters themselves");
            Assert.That(((ClassicDeconvolutionParameters)resolved.PrecursorDeconvolutionParameters).MaxAssumedChargeState,
                Is.EqualTo(12));
        });

        Assert.That(warnings, Has.Count.EqualTo(1), "A failed file-specific override must warn exactly once");
        Assert.That(warnings[0], Does.Contain(RawWithOverride), "Warning must name the raw file");
        Assert.That(warnings[0], Does.Contain(NonexistentFeaturePath), "Warning must name the selected feature path");
        Assert.That(warnings[0].ToLowerInvariant(), Does.Contain("does not exist"), "Warning must name the failure reason");
        Assert.That(warnings[0], Does.Contain(nameof(ClassicDeconvolutionParameters)),
            "Warning must name the task-wide fallback type");
    }

    [Test]
    public static void EmptyFeaturePath_WarnsAndFallsBackToTaskWideDeconvolution()
    {
        var commonParameters = ClassicCommonParameters();
        var fileSpecific = FromFileSettings(string.Empty);

        List<string> warnings = CaptureWarnings(() =>
        {
            var resolved = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                commonParameters, fileSpecific, RawWithOverride);

            Assert.That(resolved.PrecursorDeconvolutionParameters, Is.TypeOf<ClassicDeconvolutionParameters>());
        });

        Assert.That(warnings, Has.Count.EqualTo(1));
        Assert.That(warnings[0], Does.Contain(RawWithOverride));
        Assert.That(warnings[0].ToLowerInvariant(), Does.Contain("no feature file path"));
        Assert.That(warnings[0], Does.Contain(nameof(ClassicDeconvolutionParameters)));
    }

    [Test]
    public static void NonMs1ResultFile_WarnsAndFallsBackToTaskWideDeconvolution()
    {
        Assume.That(File.Exists(NonMs1FeaturePath), $"Fixture must exist at: {NonMs1FeaturePath}");

        var commonParameters = ClassicCommonParameters();
        var fileSpecific = FromFileSettings(NonMs1FeaturePath);

        List<string> warnings = CaptureWarnings(() =>
        {
            var resolved = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                commonParameters, fileSpecific, RawWithOverride);

            Assert.That(resolved.PrecursorDeconvolutionParameters, Is.TypeOf<ClassicDeconvolutionParameters>());
        });

        Assert.That(warnings, Has.Count.EqualTo(1));
        Assert.That(warnings[0], Does.Contain(RawWithOverride));
        Assert.That(warnings[0], Does.Contain(NonMs1FeaturePath));
        Assert.That(warnings[0].ToLowerInvariant(), Does.Contain("not a recognized ms1 feature file"),
            "Warning must explain that the reader result was not an MS1 feature file");
        Assert.That(warnings[0], Does.Contain(nameof(ClassicDeconvolutionParameters)));
    }

    [Test]
    public static void MalformedFeatureFile_WarnsWithParseFailureAndFallsBack()
    {
        string tempDirectory = CreateTempDirectory();
        try
        {
            string malformedPath = Path.Combine(tempDirectory, "malformed_ms1.feature");
            string header = File.ReadLines(ValidMs1FeaturePath).First();
            string[] columns = header.Split('\t');
            string row = string.Join("\t", columns.Select(column =>
                column is "Sample_ID" or "ID" ? "0"
                : column.Contains("charge_state") || column.Contains("fraction_id") ? "1"
                : "not_a_number"));
            File.WriteAllText(malformedPath, header + Environment.NewLine + row + Environment.NewLine);

            var commonParameters = ClassicCommonParameters();
            var fileSpecific = FromFileSettings(malformedPath);

            List<string> warnings = CaptureWarnings(() =>
            {
                var resolved = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                    commonParameters, fileSpecific, RawWithOverride);

                Assert.That(resolved.PrecursorDeconvolutionParameters, Is.TypeOf<ClassicDeconvolutionParameters>(),
                    "A parse failure must not escape the preflight");
            });

            Assert.That(warnings, Has.Count.EqualTo(1));
            Assert.That(warnings[0], Does.Contain(malformedPath));
            Assert.That(warnings[0].ToLowerInvariant(), Does.Contain("could not be parsed"),
                "Warning must explain the feature file could not be parsed");
            Assert.That(warnings[0], Does.Contain(nameof(ClassicDeconvolutionParameters)));
        }
        finally
        {
            DeleteTempDirectory(tempDirectory);
        }
    }

    [Test]
    public static void FileSpecificFromFile_Fallback_EmptyFeatureFile_WarnsBeforeScanProcessing()
    {
        string tempDirectory = CreateTempDirectory();
        try
        {
            string emptyPath = Path.Combine(tempDirectory, "empty_ms1.feature");
            string header = File.ReadLines(ValidMs1FeaturePath).First();
            File.WriteAllText(emptyPath, header + Environment.NewLine);

            var commonParameters = ClassicCommonParameters();
            var fileSpecific = FromFileSettings(emptyPath);

            List<string> warnings = CaptureWarnings(() =>
            {
                var resolved = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                    commonParameters, fileSpecific, RawWithOverride);

                Assert.That(resolved.PrecursorDeconvolutionParameters, Is.TypeOf<ClassicDeconvolutionParameters>(),
                    "An empty feature file must not be handed to scan processing");
            });

            Assert.That(warnings, Has.Count.EqualTo(1));
            Assert.That(warnings[0], Does.Contain(RawWithOverride));
            Assert.That(warnings[0], Does.Contain(emptyPath));
            Assert.That(warnings[0].ToLowerInvariant(), Does.Contain("no ms1 features"));
            Assert.That(warnings[0], Does.Contain(nameof(ClassicDeconvolutionParameters)));
        }
        finally
        {
            DeleteTempDirectory(tempDirectory);
        }
    }

    [Test]
    public static void ValidFromFile_IsPreloaded_AndWinsOnlyForItsRaw()
    {
        Assume.That(File.Exists(ValidMs1FeaturePath), $"Fixture must exist at: {ValidMs1FeaturePath}");

        var commonParameters = ClassicCommonParameters();
        var fileSpecific = new FileSpecificParameters
        {
            PrecursorDeconvolutionParameters = new FromFileDeconvolutionParameters(
                ValidMs1FeaturePath, minCharge: 3, maxCharge: 12, polarity: Polarity.Negative)
            {
                UseGenericScore = true
            }
        };

        List<string> warnings = CaptureWarnings(() =>
        {
            var resolvedForOverrideRaw = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                commonParameters, fileSpecific, RawWithOverride);
            var resolvedForOtherRaw = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                commonParameters, null, RawWithoutOverride);

            // The file-specific native object itself is used, so the preflight load it forced is the
            // one scan processing later consumes (no lazy surprise inside MS2 execution).
            Assert.That(resolvedForOverrideRaw.PrecursorDeconvolutionParameters,
                Is.SameAs(fileSpecific.PrecursorDeconvolutionParameters));
            var loaded = (FromFileDeconvolutionParameters)resolvedForOverrideRaw.PrecursorDeconvolutionParameters;
            Assert.That(loaded.Features, Is.Not.Empty,
                "Preflight must force the lazy feature load before scan processing");
            Assert.That(loaded.FilePath, Is.EqualTo(ValidMs1FeaturePath));
            Assert.That(loaded.MinAssumedChargeState, Is.EqualTo(3));
            Assert.That(loaded.MaxAssumedChargeState, Is.EqualTo(12));
            Assert.That(loaded.Polarity, Is.EqualTo(Polarity.Negative));
            Assert.That(loaded.UseGenericScore, Is.True);

            // A raw without an override keeps the task-wide parameters untouched.
            Assert.That(resolvedForOtherRaw.PrecursorDeconvolutionParameters,
                Is.SameAs(commonParameters.PrecursorDeconvolutionParameters));
        });

        Assert.That(warnings, Is.Empty, "A valid file-specific override must not warn");
    }

    [Test]
    public static void FileSpecificFromFile_Fallback_InvalidOverride_UsesTaskWideEmbeddedMap_ForThatRawOnly_AndOtherRawStillResolves()
    {
        Assume.That(File.Exists(ValidMs1FeaturePath), $"Fixture must exist at: {ValidMs1FeaturePath}");
        Assume.That(File.Exists(SecondValidMs1FeaturePath), $"Fixture must exist at: {SecondValidMs1FeaturePath}");

        var taskWideMap = new FeatureMappedFromFileDeconvolutionParameters(
            new SearchFeatureFileMap(new[]
            {
                new SearchFeatureFileMapEntry(RawWithOverride, ValidMs1FeaturePath),
                new SearchFeatureFileMapEntry(RawWithoutOverride, SecondValidMs1FeaturePath)
            }),
            minCharge: 1,
            maxCharge: 20)
        {
            UseGenericScore = true
        };
        var commonParameters = new CommonParameters(precursorDeconParams: taskWideMap);
        var fileSpecific = FromFileSettings(
            NonexistentFeaturePath, minCharge: 5, maxCharge: 9, polarity: Polarity.Negative);

        List<string> warnings = CaptureWarnings(() =>
        {
            var resolvedWithBadOverride = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                commonParameters, fileSpecific, RawWithOverride);
            var resolvedWithoutOverride = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                commonParameters, null, RawWithoutOverride);

            // The bad raw falls back to its own task-wide mapped route, which resolves and loads.
            Assert.That(resolvedWithBadOverride.PrecursorDeconvolutionParameters,
                Is.TypeOf<FromFileDeconvolutionParameters>());
            var mappedForBadRaw = (FromFileDeconvolutionParameters)resolvedWithBadOverride.PrecursorDeconvolutionParameters;
            Assert.That(Path.GetFullPath(mappedForBadRaw.FilePath), Is.EqualTo(Path.GetFullPath(ValidMs1FeaturePath)));
            Assert.That(mappedForBadRaw.Features, Is.Not.Empty,
                "The task-wide embedded map must resolve for the raw whose override failed");

            // The second raw is unaffected: one malformed source must not abort another raw.
            var mappedForOtherRaw = (FromFileDeconvolutionParameters)resolvedWithoutOverride.PrecursorDeconvolutionParameters;
            Assert.That(Path.GetFullPath(mappedForOtherRaw.FilePath), Is.EqualTo(Path.GetFullPath(SecondValidMs1FeaturePath)));
            Assert.That(mappedForOtherRaw.Features, Is.Not.Empty);
        });

        Assert.That(warnings, Has.Count.EqualTo(1), "Only the raw with the bad override may warn");
        Assert.That(warnings[0], Does.Contain(RawWithOverride));
        Assert.That(warnings[0], Does.Contain(NonexistentFeaturePath));
        Assert.That(warnings[0].ToLowerInvariant(), Does.Contain("does not exist"));
        Assert.That(warnings[0], Does.Contain(nameof(FeatureMappedFromFileDeconvolutionParameters)),
            "Warning must name the mapped task-wide fallback");

        Assert.That(taskWideMap.FeatureFileMap, Has.Count.EqualTo(2),
            "Resolving the fallback must not mutate the task-wide embedded map");
    }

    [Test]
    public static void InvalidOverride_DoesNotMutateCompanionToml()
    {
        string tempDirectory = CreateTempDirectory();
        try
        {
            string rawPath = Path.Combine(tempDirectory, "sample_raw.mzML");
            File.WriteAllText(rawPath, string.Empty);
            string companionTomlPath = Path.Combine(tempDirectory, "sample_raw.toml");
            File.WriteAllText(companionTomlPath,
                Toml.WriteString(FromFileSettings(NonexistentFeaturePath), MetaMorpheusTask.tomlConfig));
            byte[] tomlBefore = File.ReadAllBytes(companionTomlPath);

            // Read the companion TOML the same way RunTask does before resolving.
            var table = Toml.ReadFile(companionTomlPath, MetaMorpheusTask.tomlConfig);
            var loadedSettings = new FileSpecificParameters(table);

            List<string> warnings = CaptureWarnings(() =>
            {
                var resolved = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                    ClassicCommonParameters(), loadedSettings, rawPath);

                Assert.That(resolved.PrecursorDeconvolutionParameters, Is.TypeOf<ClassicDeconvolutionParameters>());
            });

            Assert.That(warnings, Has.Count.EqualTo(1));
            Assert.That(File.ReadAllBytes(companionTomlPath), Is.EqualTo(tomlBefore),
                "Runtime fallback must never rewrite the companion TOML");
            Assert.That(loadedSettings.PrecursorDeconvolutionParameters, Is.TypeOf<FromFileDeconvolutionParameters>(),
                "The in-memory file-specific selection must stay intact for a later retry");
            Assert.That(((FromFileDeconvolutionParameters)loadedSettings.PrecursorDeconvolutionParameters).FilePath,
                Is.EqualTo(NonexistentFeaturePath));
        }
        finally
        {
            DeleteTempDirectory(tempDirectory);
        }
    }

    [Test]
    public static void InvalidOverride_WithUnmappedTaskWideFallback_SurfacesExistingFailure()
    {
        var taskWideMap = new FeatureMappedFromFileDeconvolutionParameters(
            new SearchFeatureFileMap(new[]
            {
                new SearchFeatureFileMapEntry(RawWithoutOverride, NonexistentFeaturePath)
            }),
            minCharge: 1,
            maxCharge: 20);
        var commonParameters = new CommonParameters(precursorDeconParams: taskWideMap);
        var fileSpecific = FromFileSettings(NonexistentFeaturePath);

        List<string> warnings = CaptureWarnings(() =>
        {
            var exception = Assert.Throws<FeatureMappingException>(() =>
                MetaMorpheusTask.SetAllFileSpecificCommonParams(commonParameters, fileSpecific, RawWithOverride));

            Assert.That(exception!.Message, Does.Contain(RawWithOverride));
        });

        // The file-specific failure is still reported, and the invalid task-wide fallback surfaces its
        // own existing error instead of being swallowed by the preflight catch.
        Assert.That(warnings, Has.Count.EqualTo(1));
        Assert.That(warnings[0], Does.Contain(RawWithOverride));
        Assert.That(warnings[0], Does.Contain(nameof(FeatureMappedFromFileDeconvolutionParameters)));
    }

    [Test]
    public static void FileSpecificFromFile_Fallback_RepeatedCallsWithSameSettings_WarnOnce()
    {
        var commonParameters = ClassicCommonParameters();
        var fileSpecific = FromFileSettings(NonexistentFeaturePath);
        var results = new List<CommonParameters>();

        List<string> warnings = CaptureWarnings(() =>
        {
            for (int call = 0; call < 3; call++)
            {
                results.Add(MetaMorpheusTask.SetAllFileSpecificCommonParams(
                    commonParameters, fileSpecific, RawWithOverride));
            }
        });

        Assert.That(warnings, Has.Count.EqualTo(1),
            "Repeated resolver calls with the same file-specific settings instance must warn once");
        Assert.That(warnings[0], Does.Contain(RawWithOverride));
        Assert.That(warnings[0], Does.Contain(NonexistentFeaturePath));
        Assert.That(warnings[0].ToLowerInvariant(), Does.Contain("does not exist"));
        Assert.That(results, Has.Count.EqualTo(3));
        Assert.That(results.Select(result => result.PrecursorDeconvolutionParameters),
            Is.All.SameAs(commonParameters.PrecursorDeconvolutionParameters),
            "Every repeated call must still fall back to the task-wide precursor parameters");
    }

    [Test]
    public static void FileSpecificFromFile_Fallback_SameSettingsDifferentRaws_WarnOncePerRaw()
    {
        var commonParameters = ClassicCommonParameters();
        var fileSpecific = FromFileSettings(NonexistentFeaturePath);

        List<string> warnings = CaptureWarnings(() =>
        {
            var firstRawFirstCall = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                commonParameters, fileSpecific, RawWithOverride);
            var secondRawCall = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                commonParameters, fileSpecific, RawWithoutOverride);
            var firstRawSecondCall = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                commonParameters, fileSpecific, RawWithOverride);

            Assert.That(firstRawFirstCall.PrecursorDeconvolutionParameters,
                Is.SameAs(commonParameters.PrecursorDeconvolutionParameters));
            Assert.That(secondRawCall.PrecursorDeconvolutionParameters,
                Is.SameAs(commonParameters.PrecursorDeconvolutionParameters));
            Assert.That(firstRawSecondCall.PrecursorDeconvolutionParameters,
                Is.SameAs(commonParameters.PrecursorDeconvolutionParameters));
        });

        Assert.That(warnings, Has.Count.EqualTo(2),
            "Deduplication is per raw file, not per settings instance only");
        Assert.That(warnings.Any(warning => warning.Contains(RawWithOverride)), Is.True);
        Assert.That(warnings.Any(warning => warning.Contains(RawWithoutOverride)), Is.True);
    }

    [Test]
    public static void FileSpecificFromFile_Fallback_FreshSettingsForLaterRun_WarnAgain()
    {
        var commonParameters = ClassicCommonParameters();

        List<string> warnings = CaptureWarnings(() =>
        {
            // RunTask parses a fresh FileSpecificParameters for every task run, so the runtime
            // deduplication must not survive into the next run: a repaired file is retried and a
            // still-broken file warns again.
            MetaMorpheusTask.SetAllFileSpecificCommonParams(
                commonParameters, FromFileSettings(NonexistentFeaturePath), RawWithOverride);
            MetaMorpheusTask.SetAllFileSpecificCommonParams(
                commonParameters, FromFileSettings(NonexistentFeaturePath), RawWithOverride);
        });

        Assert.That(warnings, Has.Count.EqualTo(2),
            "A later task run with a fresh settings instance must warn again");
    }

    [Test]
    public static void FileSpecificFromFile_Fallback_ConcurrentCallsWithSameSettings_WarnOnce()
    {
        var commonParameters = ClassicCommonParameters();
        var fileSpecific = FromFileSettings(NonexistentFeaturePath);
        var warnings = new ConcurrentQueue<string>();
        EventHandler<StringEventArgs> handler = (_, e) => warnings.Enqueue(e.S);
        MetaMorpheusTask.WarnHandler += handler;
        try
        {
            var results = Enumerable.Range(0, 8)
                .AsParallel()
                .WithDegreeOfParallelism(8)
                .Select(_ => MetaMorpheusTask.SetAllFileSpecificCommonParams(
                    commonParameters, fileSpecific, RawWithOverride))
                .ToList();

            Assert.That(warnings, Has.Count.EqualTo(1),
                "Concurrent resolver calls must not warn more than once per raw/settings instance");
            Assert.That(results.Select(result => result.PrecursorDeconvolutionParameters),
                Is.All.SameAs(commonParameters.PrecursorDeconvolutionParameters));
        }
        finally
        {
            MetaMorpheusTask.WarnHandler -= handler;
        }
    }

    private static CommonParameters ClassicCommonParameters() =>
        new(precursorDeconParams: new ClassicDeconvolutionParameters(1, 12, 4, 3));

    private static FileSpecificParameters FromFileSettings(
        string featurePath, int minCharge = 2, int maxCharge = 20, Polarity polarity = Polarity.Positive) =>
        new()
        {
            PrecursorDeconvolutionParameters = new FromFileDeconvolutionParameters(
                featurePath, minCharge, maxCharge, polarity)
        };

    private static List<string> CaptureWarnings(Action action)
    {
        var warnings = new List<string>();
        EventHandler<StringEventArgs> handler = (_, e) => warnings.Add(e.S);
        MetaMorpheusTask.WarnHandler += handler;
        try
        {
            action();
        }
        finally
        {
            MetaMorpheusTask.WarnHandler -= handler;
        }

        return warnings;
    }

    private static string CreateTempDirectory()
    {
        string directory = Path.Combine(
            Path.GetTempPath(), "MetaMorpheusFromFileFallbackTests", Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(directory);
        return directory;
    }

    private static void DeleteTempDirectory(string directory)
    {
        if (Directory.Exists(directory))
        {
            Directory.Delete(directory, recursive: true);
        }
    }
}
