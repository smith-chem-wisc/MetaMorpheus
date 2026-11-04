using EngineLayer;
using EngineLayer.Deconvolution;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.IO;
using TaskLayer;

namespace Test;

/// <summary>
/// Regression coverage for the native mzLib <see cref="FromFileDeconvolutionParameters"/>
/// stored inside <see cref="FileSpecificParameters.PrecursorDeconvolutionParameters"/>.
/// These tests pin the companion-TOML contract used by the file-specific settings route:
/// the parameters serialize through <see cref="MetaMorpheusTask.tomlConfig"/> as native
/// FromFile (never as a separate top-level feature path), and a per-file override wins
/// over the task-wide embedded raw-to-feature map only for its own raw file.
/// </summary>
[TestFixture]
public static class FileSpecificFromFileTests
{
    private const string RawFileWithOverride = @"E:\data\sample1.mzML";
    private const string RawFileWithoutOverride = @"E:\data\sample2.mzML";
    private const string MappedFeaturePathForRaw1 = @"E:\features\sample1_mapped_ms1.feature";
    private const string MappedFeaturePathForRaw2 = @"E:\features\sample2_mapped_ms1.feature";

    // The file-specific path must be a real MS1 feature fixture: resolving a native FromFile override
    // forces the lazy mzLib feature load before scan processing, so a nonexistent path now falls back.
    private static string FileSpecificFeaturePath => Path.Combine(
        TestContext.CurrentContext.TestDirectory,
        "TestData", "file-specific-decon", "TaGe_SA_A549_3_snip_ms1.feature");

    [Test]
    public static void FileSpecificFromFile_TomlRoundTrip_PreservesNativeParameters()
    {
        var original = new FileSpecificParameters
        {
            PrecursorDeconvolutionParameters = new FromFileDeconvolutionParameters(
                FileSpecificFeaturePath, minCharge: 2, maxCharge: 18, polarity: Polarity.Negative)
            {
                UseGenericScore = true
            }
        };

        string toml = Toml.WriteString(original, MetaMorpheusTask.tomlConfig);
        var tomlTable = Toml.ReadString<TomlTable>(toml, MetaMorpheusTask.tomlConfig);
        var loaded = new FileSpecificParameters(tomlTable);

        Assert.That(loaded.PrecursorDeconvolutionParameters, Is.TypeOf<FromFileDeconvolutionParameters>());
        var loadedFromFile = (FromFileDeconvolutionParameters)loaded.PrecursorDeconvolutionParameters;
        Assert.That(loadedFromFile.DeconvolutionType, Is.EqualTo(DeconvolutionType.FromFile));
        Assert.That(loadedFromFile.FilePath, Is.EqualTo(FileSpecificFeaturePath));
        Assert.That(loadedFromFile.MinAssumedChargeState, Is.EqualTo(2));
        Assert.That(loadedFromFile.MaxAssumedChargeState, Is.EqualTo(18));
        Assert.That(loadedFromFile.Polarity, Is.EqualTo(Polarity.Negative));
        Assert.That(loadedFromFile.UseGenericScore, Is.True);

        // The native path is nested in the precursor table, never promoted to a top-level key.
        Assert.That(tomlTable.ContainsKey("FeatureFilePath"), Is.False,
            "Native FromFile parameters must not serialize a top-level FeatureFilePath key.");
        var precursorTable = tomlTable["PrecursorDeconvolutionParameters"].Get<TomlTable>();
        Assert.That(precursorTable.ContainsKey("FilePath"), Is.True);

        // mzLib Clone() must preserve every value that was round-tripped.
        var parameterClone = (FromFileDeconvolutionParameters)loadedFromFile.Clone();
        Assert.That(parameterClone, Is.EqualTo(loadedFromFile));
        Assert.That(parameterClone.FilePath, Is.EqualTo(FileSpecificFeaturePath));
        Assert.That(parameterClone.MinAssumedChargeState, Is.EqualTo(2));
        Assert.That(parameterClone.MaxAssumedChargeState, Is.EqualTo(18));
        Assert.That(parameterClone.Polarity, Is.EqualTo(Polarity.Negative));
        Assert.That(parameterClone.UseGenericScore, Is.True);

        // FileSpecificParameters.Clone() must keep the parsed precursor object.
        var fileSpecificClone = loaded.Clone();
        Assert.That(fileSpecificClone.PrecursorDeconvolutionParameters, Is.EqualTo(loaded.PrecursorDeconvolutionParameters));
    }

    [Test]
    public static void FileSpecificFromFile_OverridesTaskMap_ForSelectedRawOnly()
    {
        var taskWideMap = new FeatureMappedFromFileDeconvolutionParameters(
            new SearchFeatureFileMap(new[]
            {
                new SearchFeatureFileMapEntry(RawFileWithOverride, MappedFeaturePathForRaw1),
                new SearchFeatureFileMapEntry(RawFileWithoutOverride, MappedFeaturePathForRaw2)
            }),
            minCharge: 1,
            maxCharge: 20)
        {
            UseGenericScore = true
        };
        var commonParameters = new CommonParameters(precursorDeconParams: taskWideMap);
        var fileSpecific = new FileSpecificParameters
        {
            PrecursorDeconvolutionParameters = new FromFileDeconvolutionParameters(
                FileSpecificFeaturePath, minCharge: 3, maxCharge: 12, polarity: Polarity.Negative)
            {
                UseGenericScore = true
            }
        };

        var resolvedForRaw1 = MetaMorpheusTask.SetAllFileSpecificCommonParams(
            commonParameters, fileSpecific, RawFileWithOverride);
        var resolvedForRaw2 = MetaMorpheusTask.SetAllFileSpecificCommonParams(
            commonParameters, null, RawFileWithoutOverride);

        // The per-file native FromFile wins for its raw and keeps its own values.
        Assert.That(resolvedForRaw1.PrecursorDeconvolutionParameters, Is.TypeOf<FromFileDeconvolutionParameters>());
        var raw1Parameters = (FromFileDeconvolutionParameters)resolvedForRaw1.PrecursorDeconvolutionParameters;
        Assert.That(raw1Parameters.FilePath, Is.EqualTo(FileSpecificFeaturePath));
        Assert.That(raw1Parameters.MinAssumedChargeState, Is.EqualTo(3));
        Assert.That(raw1Parameters.MaxAssumedChargeState, Is.EqualTo(12));
        Assert.That(raw1Parameters.Polarity, Is.EqualTo(Polarity.Negative));
        Assert.That(raw1Parameters.UseGenericScore, Is.True);
        Assert.That(raw1Parameters.Features, Is.Not.Empty,
            "Resolving a file-specific native FromFile must force the lazy feature load");
        Assert.That(resolvedForRaw1.PrecursorDeconvolutionParameters, Is.Not.TypeOf<FeatureMappedFromFileDeconvolutionParameters>());

        // The second raw still resolves through the task-wide embedded map.
        Assert.That(resolvedForRaw2.PrecursorDeconvolutionParameters, Is.TypeOf<FromFileDeconvolutionParameters>());
        var raw2Parameters = (FromFileDeconvolutionParameters)resolvedForRaw2.PrecursorDeconvolutionParameters;
        Assert.That(raw2Parameters.FilePath, Is.EqualTo(MappedFeaturePathForRaw2));
        Assert.That(raw2Parameters.MinAssumedChargeState, Is.EqualTo(1));
        Assert.That(raw2Parameters.MaxAssumedChargeState, Is.EqualTo(20));
        Assert.That(raw2Parameters.UseGenericScore, Is.True);

        // Without the override, the first raw would have used its mapped feature file.
        var resolvedForRaw1WithoutOverride = MetaMorpheusTask.SetAllFileSpecificCommonParams(
            commonParameters, null, RawFileWithOverride);
        Assert.That(
            ((FromFileDeconvolutionParameters)resolvedForRaw1WithoutOverride.PrecursorDeconvolutionParameters).FilePath,
            Is.EqualTo(MappedFeaturePathForRaw1));

        // Resolving per-file parameters must not mutate the task-wide embedded map.
        Assert.That(taskWideMap.FeatureFileMap, Has.Count.EqualTo(2));
        Assert.That(taskWideMap.FeatureFileMap.TryGetFeaturePathForMassSpecFile(RawFileWithOverride, out var mappedPathForRaw1), Is.True);
        Assert.That(Path.GetFullPath(mappedPathForRaw1), Is.EqualTo(Path.GetFullPath(MappedFeaturePathForRaw1)));
        Assert.That(taskWideMap.FeatureFileMap.TryGetFeaturePathForMassSpecFile(RawFileWithoutOverride, out var mappedPathForRaw2), Is.True);
        Assert.That(Path.GetFullPath(mappedPathForRaw2), Is.EqualTo(Path.GetFullPath(MappedFeaturePathForRaw2)));
    }
}
