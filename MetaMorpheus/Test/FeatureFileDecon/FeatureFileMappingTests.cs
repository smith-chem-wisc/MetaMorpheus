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
    public static void TryGetFeaturePathForMassSpecFile_ReturnsFalseAndClearsPathWhenUnmapped()
    {
        var map = new SearchFeatureFileMap
        {
            Entries =
            [
                new SearchFeatureFileMapEntry(massSpecFilePath: @"E:\data\sample.mzML", featureFilePath: @"E:\features\sample.feature.tsv")
            ]
        };

        Assert.That(map.TryGetFeaturePathForMassSpecFile(@"E:\data\other.mzML", out var path), Is.False);
        Assert.That(path, Is.Empty);
    }

    [Test]
    public static void SearchFeatureFileMap_EmptyStateAndValidationCoverNullAndPopulatedEntries()
    {
        var emptyMap = new SearchFeatureFileMap();
        Assert.That(emptyMap.IsEmpty, Is.True);
        Assert.Throws<FeatureMappingException>(() => emptyMap.ValidateNotEmpty());

        var nullEntriesMap = new SearchFeatureFileMap { Entries = null };
        Assert.That(nullEntriesMap.IsEmpty, Is.True);
        Assert.Throws<FeatureMappingException>(() => nullEntriesMap.ValidateNotEmpty());

        var populatedMap = new SearchFeatureFileMap
        {
            Entries = [new SearchFeatureFileMapEntry("raw.mzML", "features.ms1.feature")]
        };
        Assert.That(populatedMap.IsEmpty, Is.False);
        Assert.DoesNotThrow(() => populatedMap.ValidateNotEmpty());
    }

    [Test]
    public static void SearchFeatureFileMap_EqualityHashAndEntryRepresentAllMappingValues()
    {
        var map = new SearchFeatureFileMap
        {
            Entries = [new SearchFeatureFileMapEntry("raw.mzML", "features.ms1.feature")]
        };
        var clone = map.Clone();
        var differentMap = new SearchFeatureFileMap
        {
            Entries = [new SearchFeatureFileMapEntry("raw.mzML", "other.ms1.feature")]
        };
        var entry = map.Entries[0];
        var entryClone = entry.Clone();
        var differentEntry = new SearchFeatureFileMapEntry("other.mzML", "features.ms1.feature");

        Assert.That(map.Equals(null), Is.False);
        Assert.That(map.Equals(map), Is.True);
        Assert.That(map.Equals(clone), Is.True);
        Assert.That(map.Equals(differentMap), Is.False);
        Assert.That(map.GetHashCode(), Is.EqualTo(clone.GetHashCode()));
        Assert.That(clone.Entries[0], Is.Not.SameAs(entry));

        Assert.That(entry.Equals(null), Is.False);
        Assert.That(entry.Equals(entry), Is.True);
        Assert.That(entry.Equals(entryClone), Is.True);
        Assert.That(((object)entry).Equals(entryClone), Is.True);
        Assert.That(((object)entry).Equals(new object()), Is.False);
        Assert.That(entry.Equals(differentEntry), Is.False);
        Assert.That(entry.GetHashCode(), Is.EqualTo(entryClone.GetHashCode()));
        Assert.That(entry.ToString(), Is.EqualTo("raw.mzML,features.ms1.feature"));
    }

    [Test]
    public static void Clone_CopiesEmbeddedMappings()
    {
        var original = CreateParameters();

        var clone = (FeatureMappedFromFileDeconvolutionParameters)original.Clone();

        Assert.That(clone, Is.EqualTo(original));
        Assert.That(clone.FeatureFileMap, Is.Not.SameAs(original.FeatureFileMap));
        Assert.That(clone.FeatureFileMap.Entries[0], Is.Not.SameAs(original.FeatureFileMap.Entries[0]));
        Assert.That(clone.MinAssumedChargeState, Is.EqualTo(original.MinAssumedChargeState));
        Assert.That(clone.MaxAssumedChargeState, Is.EqualTo(original.MaxAssumedChargeState));
        Assert.That(clone.Polarity, Is.EqualTo(original.Polarity));
        Assert.That(clone.ExpectedIsotopeSpacing, Is.EqualTo(original.ExpectedIsotopeSpacing));
        Assert.That(clone.UseGenericScore, Is.EqualTo(original.UseGenericScore));
        Assert.That(clone.ToDecoyParameters(), Is.Null);
    }

    [Test]
    public static void FeatureMappedParameters_EqualityAndHashCodeDependOnEmbeddedMapping()
    {
        var original = CreateParameters();
        var equalClone = (FeatureMappedFromFileDeconvolutionParameters)original.Clone();
        var different = new FeatureMappedFromFileDeconvolutionParameters(
            new SearchFeatureFileMap
            {
                Entries = [new SearchFeatureFileMapEntry(@"E:\data\sample.mzML", @"E:\features\other.feature.tsv")]
            },
            2,
            18,
            Polarity.Negative,
            new Averagine(),
            0.9876)
        {
            UseGenericScore = true
        };

        Assert.That(original.GetHashCode(), Is.EqualTo(equalClone.GetHashCode()));
        Assert.That(original.Equals(different), Is.False);
        Assert.That(original.Equals(new ClassicDeconvolutionParameters(1, 20, 4, 3)), Is.False);
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
    public static void TaskTomlRoundTrip_ParsesLegacyFromFileParametersWithoutEmbeddedMap()
    {
        var originalTask = new SearchTask
        {
            CommonParameters = new CommonParameters(
                precursorDeconParams: new FromFileDeconvolutionParameters(@"E:\features\legacy.feature", 2, 16, Polarity.Negative))
        };

        string toml = Toml.WriteString(originalTask, MetaMorpheusTask.tomlConfig);
        var loadedTask = Toml.ReadString<SearchTask>(toml, MetaMorpheusTask.tomlConfig);

        Assert.That(loadedTask.CommonParameters.PrecursorDeconvolutionParameters,
            Is.TypeOf<FromFileDeconvolutionParameters>());
        var loadedParameters = (FromFileDeconvolutionParameters)loadedTask.CommonParameters.PrecursorDeconvolutionParameters;
        Assert.That(loadedParameters.FilePath, Is.EqualTo(@"E:\features\legacy.feature"));
        Assert.That(loadedParameters.MinAssumedChargeState, Is.EqualTo(2));
        Assert.That(loadedParameters.MaxAssumedChargeState, Is.EqualTo(16));
        Assert.That(loadedParameters.Polarity, Is.EqualTo(Polarity.Negative));
    }

    [Test]
    public static void TaskTomlRead_RejectsMalformedSearchFeatureFileMapEntry()
    {
        var task = new SearchTask
        {
            CommonParameters = new CommonParameters(
                precursorDeconParams: new FeatureMappedFromFileDeconvolutionParameters(
                    new SearchFeatureFileMap
                    {
                        Entries = [new SearchFeatureFileMapEntry("raw.mzML", "features.ms1.feature")]
                    },
                    1,
                    20))
        };
        string toml = Toml.WriteString(task, MetaMorpheusTask.tomlConfig);
        const string serializedEntry = "raw.mzML\\tfeatures.ms1.feature";
        Assert.That(toml, Does.Contain(serializedEntry));
        toml = toml.Replace(serializedEntry, "raw.mzML", System.StringComparison.Ordinal);

        var exception = Assert.Throws<System.InvalidOperationException>(
            () => Toml.ReadString<SearchTask>(toml, MetaMorpheusTask.tomlConfig));

        Assert.That(exception!.ToString(), Does.Contain("Invalid SearchFeatureFileMapEntry"));
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
    public static void SetAllFileSpecificCommonParams_RequiresRawPathForMappedParameters()
    {
        var commonParameters = new CommonParameters(precursorDeconParams: CreateParameters());

        var exception = Assert.Throws<MetaMorpheusException>(
            () => MetaMorpheusTask.SetAllFileSpecificCommonParams(commonParameters, null, null));

        Assert.That(exception!.Message, Does.Contain("Raw file path must be provided"));
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
