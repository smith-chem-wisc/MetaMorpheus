using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.Deconvolution;
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
        var map = new SearchFeatureFileMap
        {
            Entries = inputBaseNames.Select((name, index) => new SearchFeatureFileMapEntry(massSpecFilePath: rawFiles[index],
                featureFilePath: Path.Combine(featureDirectory, name + "_ms1.feature")
            )).ToList()
        };
        var mappedParameters = new FeatureMappedFromFileDeconvolutionParameters(map, minCharge: 1, maxCharge: 20)
        {
            UseGenericScore = true
        };

        foreach (var entry in map.Entries)
        {
            Assert.That(File.Exists(entry.MassSpecFilePath), Is.True, entry.MassSpecFilePath);
            Assert.That(File.Exists(entry.FeatureFilePath), Is.True, entry.FeatureFilePath);
            var fileParameters = (FromFileDeconvolutionParameters)
                mappedParameters.ToDeconvolutionParameters(entry.MassSpecFilePath);
            Assert.That(fileParameters.Features, Is.Not.Empty, entry.FeatureFilePath);
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
        string outputDirectory = Path.Combine(temporaryRoot, "Output");
        Directory.CreateDirectory(outputDirectory);
        Directory.CreateDirectory(Path.Combine(temporaryRoot, "Task Settings"));

        try
        {
            searchTask.RunTask(
                outputDirectory,
                new List<DbForTask> { new(databasePath, false) },
                rawFiles,
                "FromFileDeconvolutionSearch");

            string psmsPath = Path.Combine(outputDirectory, "AllPSMs.psmtsv");
            Assert.That(File.Exists(psmsPath), Is.True);
            Assert.That(File.ReadAllLines(psmsPath).Length, Is.GreaterThan(1));
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
}
