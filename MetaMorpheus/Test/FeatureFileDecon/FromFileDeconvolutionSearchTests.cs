using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Deconvolution;
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

    private static string CreateRunDirectories(string runRoot)
    {
        string outputDirectory = Path.Combine(runRoot, "Output");
        Directory.CreateDirectory(outputDirectory);
        Directory.CreateDirectory(Path.Combine(runRoot, "Task Settings"));
        return outputDirectory;
    }
}
