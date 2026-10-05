using Chemistry;
using EngineLayer;
using EngineLayer.Deconvolution;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using TaskLayer;

namespace Test;

[TestFixture]
public static class FeatureFileMappingTests
{
    [Test]
    public static void TryGetFeaturePathForMassSpecFile_MatchesPathIgnoringCase()
    {
        var map = new SearchFeatureFileMap
        {
            Entries =
            [
                new SearchFeatureFileMapEntry(massSpecFilePath: @"E:\data\sample.mzML", featureFilePath: @"E:\features\sample.feature.tsv")
            ]
        };

        Assert.That(map.TryGetFeaturePathForMassSpecFile(@"e:\DATA\SAMPLE.mzml", out var path), Is.True);
        Assert.That(path, Is.EqualTo(@"E:\features\sample.feature.tsv"));
    }

    [Test]
    public static void Clone_CopiesEmbeddedMappings()
    {
        var original = CreateParameters();

        var clone = (FeatureMappedFromFileDeconvolutionParameters)original.Clone();

        Assert.That(clone, Is.EqualTo(original));
        Assert.That(clone.FeatureFileMap, Is.Not.SameAs(original.FeatureFileMap));
        Assert.That(clone.FeatureFileMap.Entries[0], Is.Not.SameAs(original.FeatureFileMap.Entries[0]));
    }

    [Test]
    public static void TaskTomlRoundTrip_PreservesMappingsAndResolvesMzLibParameters()
    {
        var originalTask = new SearchTask
        {
            CommonParameters = new CommonParameters(precursorDeconParams: CreateParameters())
        };

        string toml = Toml.WriteString(originalTask, MetaMorpheusTask.tomlConfig);
        var loadedTask = Toml.ReadString<SearchTask>(toml, MetaMorpheusTask.tomlConfig);
        var mapped = (FeatureMappedFromFileDeconvolutionParameters)loadedTask.CommonParameters.PrecursorDeconvolutionParameters;

        Assert.That(mapped.FeatureFileMap.Entries, Has.Count.EqualTo(1));
        var resolved = (FromFileDeconvolutionParameters)mapped.ToDeconvolutionParameters(@"E:\data\sample.mzML");
        Assert.That(resolved.FilePath, Is.EqualTo(@"E:\features\sample.feature.tsv"));
        Assert.That(resolved.MinAssumedChargeState, Is.EqualTo(2));
        Assert.That(resolved.MaxAssumedChargeState, Is.EqualTo(18));
        Assert.That(resolved.Polarity, Is.EqualTo(Polarity.Negative));
        Assert.That(resolved.ExpectedIsotopeSpacing, Is.EqualTo(0.9876));
        Assert.That(resolved.UseGenericScore, Is.True);
    }

    [Test]
    public static void SetAllFileSpecificCommonParams_ResolvesMappingForRawFile()
    {
        var commonParameters = new CommonParameters(precursorDeconParams: CreateParameters());

        var resolved = MetaMorpheusTask.SetAllFileSpecificCommonParams(
            commonParameters, null, @"E:\data\sample.mzML");

        Assert.That(resolved.PrecursorDeconvolutionParameters, Is.TypeOf<FromFileDeconvolutionParameters>());
        Assert.That(((FromFileDeconvolutionParameters)resolved.PrecursorDeconvolutionParameters).FilePath,
            Is.EqualTo(@"E:\features\sample.feature.tsv"));
    }

    [Test]
    public static void ToDeconvolutionParameters_ThrowsForUnmappedRawFile()
    {
        var parameters = CreateParameters();

        var exception = Assert.Throws<FeatureMappingException>(
            () => parameters.ToDeconvolutionParameters(@"E:\data\other.mzML"));

        Assert.That(exception!.Message, Does.Contain(@"E:\data\other.mzML"));
    }

    [Test]
    public static void ToDeconvolutionParameters_ThrowsForEmptyMap()
    {
        var parameters = new FeatureMappedFromFileDeconvolutionParameters();

        var exception = Assert.Throws<FeatureMappingException>(
            () => parameters.ToDeconvolutionParameters(@"E:\data\sample.mzML"));

        Assert.That(exception!.Message, Does.Contain("empty"));
    }

    private static FeatureMappedFromFileDeconvolutionParameters CreateParameters()
        => new(
            new SearchFeatureFileMap
            {
                Entries = new List<SearchFeatureFileMapEntry>
                {
                    new(massSpecFilePath: @"E:\data\sample.mzML", featureFilePath: @"E:\features\sample.feature.tsv")
                }
            },
            2,
            18,
            Polarity.Negative,
            new Averagine(),
            0.9876)
        {
            UseGenericScore = true
        };
}
