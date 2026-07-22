using Chemistry;
using EngineLayer;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using TaskLayer;
using EngineLayer.Deconvolution;
using EngineLayer.Deconvolution.FeatureFileMapping;

namespace Test;

[TestFixture]
public static class FeatureFileMappingTests
{
    [Test]
    public static void BuildSearchMap_IncludesOnlyRequestedConditionAndSearchFiles()
    {
        var map = new FeatureFileMap
        {
            FilePath = @"E:\maps\feature-map-a.toml"
        };
        map.AddOrReplaceCondition(new FeatureFileMapCondition
        {
            ConditionKey = "flash",
            DisplayName = "FlashDeconv"
        });
        map.AddOrReplaceCondition(new FeatureFileMapCondition
        {
            ConditionKey = "dino",
            DisplayName = "Dinosaur"
        });

        var sample1 = FeatureSpectraEntry.Create(@"E:\data\sample1.mzML");
        sample1.AddOrReplaceConditionFile(FeatureFileConditionEntry.Create("flash", @"E:\features\flash\sample1_ms1.feature"));
        sample1.AddOrReplaceConditionFile(FeatureFileConditionEntry.Create("dino", @"E:\features\dino\sample1.feature.tsv"));
        map.AddOrReplaceSpectraEntry(sample1);

        var sample2 = FeatureSpectraEntry.Create(@"E:\data\sample2.mzML");
        sample2.AddOrReplaceConditionFile(FeatureFileConditionEntry.Create("flash", @"E:\features\flash\sample2_ms1.feature"));
        map.AddOrReplaceSpectraEntry(sample2);

        var searchMap = map.BuildSearchMap(new[] { @"E:\data\sample1.mzML" }, "flash");

        Assert.That(searchMap.SourceMapPath, Is.EqualTo(map.FilePath));
        Assert.That(searchMap.SelectedConditionKey, Is.EqualTo("flash"));
        Assert.That(searchMap.SelectedConditionDisplayName, Is.EqualTo("FlashDeconv"));
        Assert.That(searchMap.Entries, Has.Count.EqualTo(1));
        Assert.That(searchMap.Entries[0].MassSpecFilePath, Is.EqualTo(@"E:\data\sample1.mzML"));
        Assert.That(searchMap.Entries[0].MassSpecFileName, Is.EqualTo("sample1.mzML"));
        Assert.That(searchMap.Entries[0].FeatureFilePath, Is.EqualTo(@"E:\features\flash\sample1_ms1.feature"));
        Assert.That(searchMap.Entries[0].FeatureFileName, Is.EqualTo("sample1_ms1.feature"));
    }

    [Test]
    public static void BuildSearchMap_ThrowsInformativeExceptionWhenConditionMissing()
    {
        var map = new FeatureFileMap();

        var exception = Assert.Throws<FeatureMappingException>(() => map.BuildSearchMap(new[] { @"E:\data\sample1.mzML" }, "flash"));

        Assert.That(exception!.Message, Does.Contain("Condition 'flash' not found in map."));
    }

    [Test]
    public static void BuildSearchMap_ThrowsInformativeExceptionWhenRawFileMissingConditionMapping()
    {
        var map = new FeatureFileMap();
        map.AddOrReplaceCondition(new FeatureFileMapCondition { ConditionKey = "flash", DisplayName = "FlashDeconv" });

        var sample1 = FeatureSpectraEntry.Create(@"E:\data\sample1.mzML");
        sample1.AddOrReplaceConditionFile(FeatureFileConditionEntry.Create("flash", @"E:\features\flash\sample1_ms1.feature"));
        map.AddOrReplaceSpectraEntry(sample1);

        var exception = Assert.Throws<FeatureMappingException>(() => map.BuildSearchMap(
            new[] { @"E:\data\sample1.mzML", @"E:\data\sample2.mzML" },
            "flash"));

        Assert.That(exception!.Message, Does.Contain(@"E:\data\sample2.mzML"));
        Assert.That(exception.Message, Does.Contain("condition 'flash'"));
    }

    [Test]
    public static void ToDeconvolutionParameters_CopiesInheritedDeconvolutionSettings()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);

        var searchMap = new SearchFeatureFileMap
        {
            Entries = new List<SearchFeatureFileMapEntry>
            {
                new()
                {
                    MassSpecFilePath = Path.Combine(tempDir, "sample1.mzML"),
                    MassSpecFileName = "sample1.mzML",
                    FeatureFilePath = Path.Combine(tempDir, "sample1.feature.tsv"),
                    FeatureFileName = "sample1.feature.tsv"
                }
            }
        };

        File.WriteAllText(searchMap.Entries[0].FeatureFilePath, string.Empty);

        var mappedParameters = new FeatureMappedFromFileDeconvolutionParameters(searchMap, 2, 18, Polarity.Negative, new Averagine(), 0.9876)
        {
            UseGenericScore = true
        };

        try
        {
            var resolved = (FromFileDeconvolutionParameters)mappedParameters.ToDeconvolutionParameters(searchMap.Entries[0].MassSpecFilePath);
            Assert.That(resolved.MinAssumedChargeState, Is.EqualTo(2));
            Assert.That(resolved.MaxAssumedChargeState, Is.EqualTo(18));
            Assert.That(resolved.Polarity, Is.EqualTo(Polarity.Negative));
            Assert.That(resolved.ExpectedIsotopeSpacing, Is.EqualTo(0.9876));
            Assert.That(resolved.UseGenericScore, Is.True);
            Assert.That(resolved.AverageResidueModel, Is.TypeOf<Averagine>());
        }
        finally
        {
            if (Directory.Exists(tempDir))
            {
                Directory.Delete(tempDir, true);
            }
        }
    }

    [Test]
    public static void SearchFeatureFileMapClone_DeepClonesEntries()
    {
        var map = new SearchFeatureFileMap
        {
            SourceMapPath = @"E:\maps\feature-map-a.toml",
            SelectedConditionKey = "flash",
            SelectedConditionDisplayName = "FlashDeconv",
            Entries = new List<SearchFeatureFileMapEntry>
            {
                new()
                {
                    MassSpecFilePath = @"E:\data\sample1.mzML",
                    MassSpecFileName = "sample1.mzML",
                    FeatureFilePath = @"E:\features\sample1_ms1.feature",
                    FeatureFileName = "sample1_ms1.feature"
                }
            }
        };

        var clone = map.Clone();
        clone.Entries[0].FeatureFilePath = @"E:\features\changed.feature";

        Assert.That(map.Entries[0].FeatureFilePath, Is.EqualTo(@"E:\features\sample1_ms1.feature"));
        Assert.That(clone.Entries[0].FeatureFilePath, Is.EqualTo(@"E:\features\changed.feature"));
    }

    [Test]
    public static void SearchFeatureFileMap_Equality_IsValueBased()
    {
        var left = new SearchFeatureFileMap
        {
            SourceMapPath = @"E:\maps\feature-map-a.toml",
            SelectedConditionKey = "flash",
            SelectedConditionDisplayName = "FlashDeconv",
            Entries = new List<SearchFeatureFileMapEntry>
            {
                new()
                {
                    MassSpecFilePath = @"E:\data\sample1.mzML",
                    MassSpecFileName = "sample1.mzML",
                    FeatureFilePath = @"E:\features\sample1_ms1.feature",
                    FeatureFileName = "sample1_ms1.feature"
                }
            }
        };

        var right = left.Clone();

        Assert.That(left, Is.EqualTo(right));
        Assert.That(left.GetHashCode(), Is.EqualTo(right.GetHashCode()));
    }

    [Test]
    public static void FeatureMappedFromFileDeconvolutionParameters_Equality_UsesSearchMapValueEquality()
    {
        var searchMap = new SearchFeatureFileMap
        {
            SourceMapPath = @"E:\maps\feature-map-a.toml",
            SelectedConditionKey = "flash",
            SelectedConditionDisplayName = "FlashDeconv",
            Entries = new List<SearchFeatureFileMapEntry>
            {
                new()
                {
                    MassSpecFilePath = @"E:\data\sample1.mzML",
                    MassSpecFileName = "sample1.mzML",
                    FeatureFilePath = @"E:\features\sample1_ms1.feature",
                    FeatureFileName = "sample1_ms1.feature"
                }
            }
        };

        var left = new FeatureMappedFromFileDeconvolutionParameters(searchMap, 2, 18, Polarity.Positive, new Averagine(), 0.9876)
        {
            UseGenericScore = true
        };
        var right = new FeatureMappedFromFileDeconvolutionParameters(searchMap.Clone(), 2, 18, Polarity.Positive, new Averagine(), 0.9876)
        {
            UseGenericScore = true
        };

        Assert.That(left, Is.EqualTo(right));
        Assert.That(left.GetHashCode(), Is.EqualTo(right.GetHashCode()));
    }

    [Test]
    public static void FeatureMappedFromFileDeconvolutionParameters_TomlRoundTrip_PreservesSearchMap()
    {
        var original = new FeatureMappedFromFileDeconvolutionParameters(
            new SearchFeatureFileMap
            {
                SourceMapPath = @"E:\maps\feature-map-a.toml",
                SelectedConditionKey = "flash",
                SelectedConditionDisplayName = "FlashDeconv",
                Entries = new List<SearchFeatureFileMapEntry>
                {
                    new()
                    {
                        MassSpecFilePath = @"E:\data\sample1.mzML",
                        MassSpecFileName = "sample1.mzML",
                        FeatureFilePath = @"E:\features\sample1_ms1.feature",
                        FeatureFileName = "sample1_ms1.feature"
                    },
                    new()
                    {
                        MassSpecFilePath = @"E:\data\sample2.mzML",
                        MassSpecFileName = "sample2.mzML",
                        FeatureFilePath = @"E:\features\sample2_ms1.feature",
                        FeatureFileName = "sample2_ms1.feature"
                    }
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

        var searchTask = new SearchTask
        {
            CommonParameters = new CommonParameters(precursorDeconParams: original)
        };

        string toml = Toml.WriteString(searchTask, MetaMorpheusTask.tomlConfig);
        var roundTripped = Toml.ReadString<SearchTask>(toml, MetaMorpheusTask.tomlConfig);
        var parsed = roundTripped.CommonParameters.PrecursorDeconvolutionParameters as FeatureMappedFromFileDeconvolutionParameters;

        Assert.That(parsed, Is.Not.Null);
        Assert.That(parsed, Is.EqualTo(original));
        Assert.That(parsed!.FeatureFileMap, Is.EqualTo(original.FeatureFileMap));
    }

    [Test]
    public static void SetAllFileSpecificCommonParams_MaterializesMappedPrecursorPerRawFile()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);

        try
        {
            string rawFilePath = Path.Combine(tempDir, "sample1.mzML");
            string featureFilePath = Path.Combine(tempDir, "sample1.feature.tsv");
            File.WriteAllText(featureFilePath, string.Empty);

            var mappedParameters = new FeatureMappedFromFileDeconvolutionParameters(
                new SearchFeatureFileMap
                {
                    Entries = new List<SearchFeatureFileMapEntry>
                    {
                        new()
                        {
                            MassSpecFilePath = rawFilePath,
                            MassSpecFileName = "sample1.mzML",
                            FeatureFilePath = featureFilePath,
                            FeatureFileName = "sample1.feature.tsv"
                        }
                    }
                },
                3,
                17,
                Polarity.Positive,
                new Averagine(),
                0.9988)
            {
                UseGenericScore = true
            };

            var commonParameters = new CommonParameters(
                precursorDeconParams: mappedParameters,
                productDeconParams: new ClassicDeconvolutionParameters(1, 10, 4, 3));

            var resolved = MetaMorpheusTask.SetAllFileSpecificCommonParams(commonParameters, null, rawFilePath);

            Assert.That(resolved.PrecursorDeconvolutionParameters, Is.TypeOf<FromFileDeconvolutionParameters>());
            var precursor = (FromFileDeconvolutionParameters)resolved.PrecursorDeconvolutionParameters;
            Assert.That(precursor.MinAssumedChargeState, Is.EqualTo(3));
            Assert.That(precursor.MaxAssumedChargeState, Is.EqualTo(17));
            Assert.That(precursor.ExpectedIsotopeSpacing, Is.EqualTo(0.9988));
            Assert.That(precursor.UseGenericScore, Is.True);
            Assert.That(resolved.ProductDeconvolutionParameters, Is.TypeOf<ClassicDeconvolutionParameters>());
        }
        finally
        {
            if (Directory.Exists(tempDir))
            {
                Directory.Delete(tempDir, true);
            }
        }
    }

    [Test]
    public static void ToDeconvolutionParameters_ThrowsOnEmptyMap()
    {
        var emptyParams = new FeatureMappedFromFileDeconvolutionParameters();

        var ex = Assert.Throws<FeatureMappingException>(() =>
            emptyParams.ToDeconvolutionParameters(@"E:\data\missing.mzML"));

        Assert.That(ex!.Message, Does.Contain("empty"));
        Assert.That(ex.Message, Does.Contain("no feature file entries"));
    }

    [Test]
    public static void SearchFeatureFileMap_Equality_IgnoresSourceMapPath()
    {
        var entries = new List<SearchFeatureFileMapEntry>
        {
            new()
            {
                MassSpecFilePath = @"E:\data\a.mzML",
                MassSpecFileName = "a.mzML",
                FeatureFilePath = @"E:\features\a.feature",
                FeatureFileName = "a.feature"
            }
        };

        var mapA = new SearchFeatureFileMap
        {
            SourceMapPath = @"E:\maps\original.toml",
            SelectedConditionKey = "flash",
            SelectedConditionDisplayName = "FlashDeconv",
            Entries = entries
        };

        var mapB = new SearchFeatureFileMap
        {
            SourceMapPath = @"E:\maps\different-path.toml",
            SelectedConditionKey = "flash",
            SelectedConditionDisplayName = "FlashDeconv",
            Entries = entries
        };

        // Different SourceMapPath must not break equality
        Assert.That(mapA, Is.EqualTo(mapB));
        Assert.That(mapA.GetHashCode(), Is.EqualTo(mapB.GetHashCode()));
    }

    [Test]
    public static void SearchFeatureFileMap_IsEmpty_TrueWhenNoEntries()
    {
        var map = new SearchFeatureFileMap();
        Assert.That(map.IsEmpty, Is.True);
    }

    [Test]
    public static void SearchFeatureFileMap_IsEmpty_FalseWhenHasEntries()
    {
        var map = new SearchFeatureFileMap
        {
            Entries = new List<SearchFeatureFileMapEntry>
            {
                new()
                {
                    MassSpecFilePath = @"E:\data\a.mzML",
                    FeatureFilePath = @"E:\features\a.feature"
                }
            }
        };
        Assert.That(map.IsEmpty, Is.False);
    }

    [Test]
    public static void ValidateNotEmpty_ThrowsOnEmptyMap()
    {
        var map = new SearchFeatureFileMap();

        var ex = Assert.Throws<FeatureMappingException>(() => map.ValidateNotEmpty());

        Assert.That(ex!.Message, Does.Contain("empty"));
    }

    [Test]
    public static void ValidateNotEmpty_PassesOnNonEmptyMap()
    {
        var map = new SearchFeatureFileMap
        {
            Entries = new List<SearchFeatureFileMapEntry>
            {
                new()
                {
                    MassSpecFilePath = @"E:\data\a.mzML",
                    FeatureFilePath = @"E:\features\a.feature"
                }
            }
        };

        Assert.DoesNotThrow(() => map.ValidateNotEmpty());
    }

    [Test]
    public static void SourceMapPath_IsPreservedInClone()
    {
        var map = new SearchFeatureFileMap
        {
            SourceMapPath = @"E:\maps\provenance.toml",
            SelectedConditionKey = "dino",
            Entries = new List<SearchFeatureFileMapEntry>
            {
                new()
                {
                    MassSpecFilePath = @"E:\data\x.mzML",
                    FeatureFilePath = @"E:\features\x.feature"
                }
            }
        };

        var clone = map.Clone();

        // SourceMapPath must survive clone (audit trail preservation)
        Assert.That(clone.SourceMapPath, Is.EqualTo(@"E:\maps\provenance.toml"));
    }

    [Test]
    public static void SetAllFileSpecificCommonParams_ThrowsOnEmptyMappedDeconvParams()
    {
        var emptyMapped = new FeatureMappedFromFileDeconvolutionParameters();
        var commonParams = new CommonParameters(precursorDeconParams: emptyMapped);

        var ex = Assert.Throws<FeatureMappingException>(() =>
            MetaMorpheusTask.SetAllFileSpecificCommonParams(commonParams, null, @"E:\data\sample.mzML"));

        Assert.That(ex!.Message, Does.Contain("empty"));
    }

    #region FeatureFileMapStore tests

    [Test]
    public static void Store_Load_ReturnsEmptyMapWhenFileMissing()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            string missingPath = Path.Combine(tempDir, "does-not-exist.toml");
            var map = FeatureFileMapStore.Load(missingPath);

            Assert.That(map, Is.Not.Null);
            Assert.That(map.FilePath, Is.EqualTo(missingPath));
            Assert.That(map.SpectraFiles, Is.Empty);
            Assert.That(map.Conditions, Is.Empty);
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public static void Store_SaveAndLoad_RoundTripsContent()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            string storePath = Path.Combine(tempDir, "feature-maps.toml");

            var original = new FeatureFileMap();
            original.AddOrReplaceCondition(new FeatureFileMapCondition
            {
                ConditionKey = "flash",
                DisplayName = "FlashDeconv"
            });
            var entry = FeatureSpectraEntry.Create(@"E:\data\sample.mzML");
            entry.AddOrReplaceConditionFile(FeatureFileConditionEntry.Create("flash", @"E:\features\sample.feature"));
            original.AddOrReplaceSpectraEntry(entry);

            FeatureFileMapStore.Save(original, storePath);

            Assert.That(File.Exists(storePath), Is.True);

            var loaded = FeatureFileMapStore.Load(storePath);

            Assert.That(loaded.FilePath, Is.EqualTo(storePath));
            Assert.That(loaded.Conditions, Has.Count.EqualTo(1));
            Assert.That(loaded.Conditions[0].ConditionKey, Is.EqualTo("flash"));
            Assert.That(loaded.SpectraFiles, Has.Count.EqualTo(1));
            Assert.That(loaded.SpectraFiles[0].MassSpecFilePath, Is.EqualTo(@"E:\data\sample.mzML"));

            // Verify condition file came through
            loaded.TryGetSpectraEntry(@"E:\data\sample.mzML", out var spectraEntry);
            Assert.That(spectraEntry, Is.Not.Null);
            spectraEntry.TryGetConditionFile("flash", out var conditionFile);
            Assert.That(conditionFile, Is.Not.Null);
            Assert.That(conditionFile.FeatureFilePath, Is.EqualTo(@"E:\features\sample.feature"));
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public static void Store_Load_ThrowsOnInvalidToml()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            string storePath = Path.Combine(tempDir, "corrupted.toml");
            File.WriteAllText(storePath, "<<<NOT VALID TOML>>>");

            var ex = Assert.Throws<FeatureMappingException>(() => FeatureFileMapStore.Load(storePath));

            Assert.That(ex!.Message, Does.Contain("Failed to parse feature-map TOML"));
            Assert.That(ex.Message, Does.Contain(storePath));
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public static void Store_Save_AtomicReplacePreservesContentOnRewrite()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            string storePath = Path.Combine(tempDir, "feature-maps.toml");

            // First save
            var first = new FeatureFileMap();
            first.AddOrReplaceCondition(new FeatureFileMapCondition { ConditionKey = "v1" });
            first.AddOrReplaceCondition(new FeatureFileMapCondition { ConditionKey = "v2" });
            FeatureFileMapStore.Save(first, storePath);

            // Second save (replaces existing via File.Replace)
            var second = new FeatureFileMap();
            second.AddOrReplaceCondition(new FeatureFileMapCondition { ConditionKey = "v3" });
            FeatureFileMapStore.Save(second, storePath);

            // Load back — should have "v3" only
            var loaded = FeatureFileMapStore.Load(storePath);
            Assert.That(loaded.Conditions, Has.Count.EqualTo(1));
            Assert.That(loaded.Conditions[0].ConditionKey, Is.EqualTo("v3"));

            // No orphaned .tmp.* files should remain
            var tmpFiles = Directory.GetFiles(tempDir, "*.tmp.*");
            Assert.That(tmpFiles, Is.Empty);
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public static void Store_Load_MissingFileSetsLastModifiedUtc()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            string missingPath = Path.Combine(tempDir, "missing.toml");
            var before = DateTime.UtcNow.AddSeconds(-1);
            var map = FeatureFileMapStore.Load(missingPath);
            var after = DateTime.UtcNow.AddSeconds(1);

            Assert.That(map.LastModifiedUtc, Is.GreaterThan(before));
            Assert.That(map.LastModifiedUtc, Is.LessThan(after));
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public static void Store_Save_ThrowsOnNullMap()
    {
        var ex = Assert.Throws<ArgumentNullException>(() => FeatureFileMapStore.Save(null));
        Assert.That(ex!.ParamName, Is.EqualTo("map"));
    }

    [Test]
    public static void Store_Load_DefaultPathDelegatesToGlobalVariables()
    {
        // The default store path must match GlobalVariables.FeatureMapsFilePath
        // This test verifies the store delegates to the singleton when no path is given
        var map = FeatureFileMapStore.Load();
        Assert.That(map, Is.Not.Null);
        Assert.That(map.FilePath, Is.EqualTo(EngineLayer.GlobalVariables.FeatureMapsFilePath));
    }

    #endregion

    #region Path normalization tests

    [Test]
    public static void Store_CreateEntry_NormalizesForwardSlashToBackslash()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            // Use forward-slash path (non-native on Windows)
            string forwardSlashPath = tempDir.Replace("\\", "/") + "/sample.mzML";
            string expectedNormalized = Path.GetFullPath(forwardSlashPath);

            var entry = FeatureSpectraEntry.Create(forwardSlashPath);

            Assert.That(entry.MassSpecFilePath, Is.EqualTo(expectedNormalized));
            Assert.That(entry.MassSpecFileName, Is.EqualTo("sample.mzML"));
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public static void Store_CreateConditionEntry_NormalizesFeatureFilePath()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            string forwardSlashPath = tempDir.Replace("\\", "/") + "/features.tsv";
            string expectedNormalized = Path.GetFullPath(forwardSlashPath);

            var cf = FeatureFileConditionEntry.Create("flash", forwardSlashPath);

            Assert.That(cf.FeatureFilePath, Is.EqualTo(expectedNormalized));
            Assert.That(cf.FeatureFileName, Is.EqualTo("features.tsv"));
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public static void Store_MatchesRawFile_EquivalentPaths_DifferentCase()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            string canonical = Path.Combine(tempDir, "Sample.mzML");
            string differentCase = Path.Combine(tempDir, "SAMPLE.mzML");

            var entry = FeatureSpectraEntry.Create(canonical);

            Assert.That(entry.MatchesRawFile(differentCase), Is.True,
                "MatchesRawFile should return true for paths differing only in casing");
            Assert.That(entry.MatchesRawFile(canonical), Is.True,
                "MatchesRawFile should return true for the exact same path");
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public static void Store_MatchesRawFile_EquivalentPaths_DifferentSeparator()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            string backslashPath = Path.Combine(tempDir, "sample.mzML");
            string forwardSlashPath = tempDir.Replace("\\", "/") + "/sample.mzML";

            var entry = FeatureSpectraEntry.Create(backslashPath);

            // Both paths resolve to the same canonical form
            Assert.That(entry.MatchesRawFile(forwardSlashPath), Is.True,
                "MatchesRawFile should normalize forward-slash paths to match backslash-stored paths");
            Assert.That(entry.MatchesRawFile(backslashPath), Is.True,
                "MatchesRawFile should match the same backslash path");
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public static void Store_Load_NormalizesPathsFromDeserializedToml()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            string storePath = Path.Combine(tempDir, "feature-maps.toml");

            // Build a map with forward-slash paths, save it
            string nonNormalizedRaw = tempDir.Replace("\\", "/") + "/sample.mzML";
            string nonNormalizedFeat = tempDir.Replace("\\", "/") + "/sample.feature";

            var map = new FeatureFileMap();
            var entry = FeatureSpectraEntry.Create(nonNormalizedRaw);
            entry.AddOrReplaceConditionFile(FeatureFileConditionEntry.Create("flash", nonNormalizedFeat));
            map.AddOrReplaceSpectraEntry(entry);
            FeatureFileMapStore.Save(map, storePath);

            // Load back — paths must be normalized
            var loaded = FeatureFileMapStore.Load(storePath);

            Assert.That(loaded.SpectraFiles[0].MassSpecFilePath,
                Is.EqualTo(Path.GetFullPath(nonNormalizedRaw)));
            Assert.That(loaded.SpectraFiles[0].MassSpecFileName, Is.EqualTo("sample.mzML"));

            loaded.SpectraFiles[0].TryGetConditionFile("flash", out var cf);
            Assert.That(cf, Is.Not.Null);
            Assert.That(cf.FeatureFilePath,
                Is.EqualTo(Path.GetFullPath(nonNormalizedFeat)));
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public static void Store_MatchesRawFile_NullOrEmptyInput()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            var emptyEntry = new FeatureSpectraEntry(); // default MassSpecFilePath = ""
            Assert.That(emptyEntry.MatchesRawFile(null), Is.True,
                "null input should match empty stored path");
            Assert.That(emptyEntry.MatchesRawFile(""), Is.True,
                "empty input should match empty stored path");

            var nonEmptyEntry = FeatureSpectraEntry.Create(
                Path.Combine(tempDir, "real.mzML"));
            Assert.That(nonEmptyEntry.MatchesRawFile(null), Is.False,
                "null input should not match non-empty stored path");
            Assert.That(nonEmptyEntry.MatchesRawFile(""), Is.False,
                "empty input should not match non-empty stored path");
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public static void GUI_and_CLI_Toml_Reproducibility_Test()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            string rawFilePath = Path.Combine(tempDir, "sample1.mzML");
            string featureFilePath = Path.Combine(tempDir, "sample1.feature.tsv");
            File.WriteAllText(featureFilePath, "mz\tmostAbundantMz\tcharge\trtStart\trtApex\trtEnd\tfwhm\tnIsotopes\tnScans\taveragineCorr\tmass\tmassCalib\tintensityApex\tintensitySum");

            // 1. GUI creates parameters and exports task
            var guiMap = new SearchFeatureFileMap
            {
                SourceMapPath = @"E:\maps\feature-map-a.toml",
                SelectedConditionKey = "flash",
                SelectedConditionDisplayName = "FlashDeconv",
                Entries = new List<SearchFeatureFileMapEntry>
                {
                    new()
                    {
                        MassSpecFilePath = rawFilePath,
                        MassSpecFileName = "sample1.mzML",
                        FeatureFilePath = featureFilePath,
                        FeatureFileName = "sample1.feature.tsv"
                    }
                }
            };
            
            var originalParams = new FeatureMappedFromFileDeconvolutionParameters(guiMap, 2, 18, Polarity.Negative, new Averagine(), 0.9876) { UseGenericScore = true };
            var originalTask = new SearchTask { CommonParameters = new CommonParameters(precursorDeconParams: originalParams) };
            
            string taskTomlPath = Path.Combine(tempDir, "task.toml");
            Toml.WriteFile(originalTask, taskTomlPath, MetaMorpheusTask.tomlConfig);

            // 2. Simulate CLI (or fresh GUI reopen) - load TOML and materialize for raw file
            // No global store is available (GlobalVariables is not even pointing here)
            var loadedTask = Toml.ReadFile<SearchTask>(taskTomlPath, MetaMorpheusTask.tomlConfig);
            var loadedCommonParams = loadedTask.CommonParameters;
            
            // Prove global store access is not needed by directly invoking materialization
            var fileSpecificParams = MetaMorpheusTask.SetAllFileSpecificCommonParams(loadedCommonParams, null, rawFilePath);
            
            Assert.That(fileSpecificParams.PrecursorDeconvolutionParameters, Is.TypeOf<FromFileDeconvolutionParameters>());
            var precursor = (FromFileDeconvolutionParameters)fileSpecificParams.PrecursorDeconvolutionParameters;
            Assert.That(precursor.MinAssumedChargeState, Is.EqualTo(2));
            Assert.That(precursor.MaxAssumedChargeState, Is.EqualTo(18));
        }
        finally
        {
            if (Directory.Exists(tempDir))
            {
                Directory.Delete(tempDir, true);
            }
        }
    }

    #endregion

    #region T10 Reproducibility Tests

    [Test]
    public static void EmbeddedMap_AllEntriesMaterializePerFile_AfterTomlRoundTrip()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            // Create feature files for two spectra files (full header required by format reader)
            string header = "mz\tmostAbundantMz\tcharge\trtStart\trtApex\trtEnd\tfwhm\tnIsotopes\tnScans\taveragineCorr\tmass\tmassCalib\tintensityApex\tintensitySum";
            string feat1Path = Path.Combine(tempDir, "sample1.feature.tsv");
            File.WriteAllText(feat1Path, header);
            string feat2Path = Path.Combine(tempDir, "sample2.feature.tsv");
            File.WriteAllText(feat2Path, header);

            string raw1Path = Path.Combine(tempDir, "sample1.mzML");
            string raw2Path = Path.Combine(tempDir, "sample2.mzML");

            // Build embedded map with TWO entries
            var map = new SearchFeatureFileMap
            {
                SourceMapPath = @"E:\maps\nonexistent-store.toml",
                SelectedConditionKey = "cond",
                SelectedConditionDisplayName = "Condition",
                Entries = new List<SearchFeatureFileMapEntry>
                {
                    new()
                    {
                        MassSpecFilePath = raw1Path,
                        MassSpecFileName = "sample1.mzML",
                        FeatureFilePath = feat1Path,
                        FeatureFileName = "sample1.feature.tsv"
                    },
                    new()
                    {
                        MassSpecFilePath = raw2Path,
                        MassSpecFileName = "sample2.mzML",
                        FeatureFilePath = feat2Path,
                        FeatureFileName = "sample2.feature.tsv"
                    }
                }
            };

            var originalParams = new FeatureMappedFromFileDeconvolutionParameters(map, 2, 18, Polarity.Positive, new Averagine(), 0.9876) { UseGenericScore = true };
            var originalTask = new SearchTask { CommonParameters = new CommonParameters(precursorDeconParams: originalParams) };

            // TOML round-trip: save → load
            string tomlPath = Path.Combine(tempDir, "task.toml");
            Toml.WriteFile(originalTask, tomlPath, MetaMorpheusTask.tomlConfig);
            var loadedTask = Toml.ReadFile<SearchTask>(tomlPath, MetaMorpheusTask.tomlConfig);
            var loadedMapped = loadedTask.CommonParameters.PrecursorDeconvolutionParameters as FeatureMappedFromFileDeconvolutionParameters;
            Assert.That(loadedMapped, Is.Not.Null);

            // Prove per-file materialization works for EVERY entry in the embedded map
            var resolved1 = loadedMapped.ToDeconvolutionParameters(raw1Path);
            Assert.That(resolved1, Is.TypeOf<FromFileDeconvolutionParameters>());
            Assert.That(((FromFileDeconvolutionParameters)resolved1).MinAssumedChargeState, Is.EqualTo(2));
            Assert.That(((FromFileDeconvolutionParameters)resolved1).MaxAssumedChargeState, Is.EqualTo(18));

            var resolved2 = loadedMapped.ToDeconvolutionParameters(raw2Path);
            Assert.That(resolved2, Is.TypeOf<FromFileDeconvolutionParameters>());
            Assert.That(((FromFileDeconvolutionParameters)resolved2).MinAssumedChargeState, Is.EqualTo(2));
            Assert.That(((FromFileDeconvolutionParameters)resolved2).MaxAssumedChargeState, Is.EqualTo(18));

            // The two resolved objects are distinct instances (not cached/shared)
            Assert.That(resolved1, Is.Not.SameAs(resolved2));

            // The full embedded map survives round-trip (value equality)
            Assert.That(loadedMapped.FeatureFileMap, Is.EqualTo(originalParams.FeatureFileMap));
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    [Test]
    public static void EmbeddedMapIsAuthoritative_SourceMapPathIsAuditOnly()
    {
        string tempDir = Path.Combine(TestContext.CurrentContext.TestDirectory, Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(tempDir);
        try
        {
            string rawPath = Path.Combine(tempDir, "sample.mzML");
            string featPath = Path.Combine(tempDir, "sample.feature.tsv");
            string header = "mz\tmostAbundantMz\tcharge\trtStart\trtApex\trtEnd\tfwhm\tnIsotopes\tnScans\taveragineCorr\tmass\tmassCalib\tintensityApex\tintensitySum";
            File.WriteAllText(featPath, header);

            // SourceMapPath points to a non-existent file — the global store is absent
            var map = new SearchFeatureFileMap
            {
                SourceMapPath = Path.Combine(tempDir, "no-such-store.toml"),
                SelectedConditionKey = "cond",
                SelectedConditionDisplayName = "Condition",
                Entries = new List<SearchFeatureFileMapEntry>
                {
                    new()
                    {
                        MassSpecFilePath = rawPath,
                        MassSpecFileName = "sample.mzML",
                        FeatureFilePath = featPath,
                        FeatureFileName = "sample.feature.tsv"
                    }
                }
            };

            Assert.That(File.Exists(map.SourceMapPath), Is.False,
                "Precondition: SourceMapPath must not exist to prove no store dependency");

            var originalParams = new FeatureMappedFromFileDeconvolutionParameters(map, 3, 20, Polarity.Negative, new Averagine(), 1.0);
            var originalTask = new SearchTask { CommonParameters = new CommonParameters(precursorDeconParams: originalParams) };

            // TOML round-trip
            string tomlPath = Path.Combine(tempDir, "task.toml");
            Toml.WriteFile(originalTask, tomlPath, MetaMorpheusTask.tomlConfig);
            var loadedTask = Toml.ReadFile<SearchTask>(tomlPath, MetaMorpheusTask.tomlConfig);
            var loadedMapped = loadedTask.CommonParameters.PrecursorDeconvolutionParameters as FeatureMappedFromFileDeconvolutionParameters;
            Assert.That(loadedMapped, Is.Not.Null);

            // SourceMapPath is preserved through serialization (audit trail works)
            Assert.That(loadedMapped.FeatureFileMap.SourceMapPath, Is.EqualTo(map.SourceMapPath),
                "SourceMapPath audit trail must survive TOML round-trip");

            // Materialization via ToDeconvolutionParameters uses embedded entries, NOT the store
            var resolved = loadedMapped.ToDeconvolutionParameters(rawPath);
            Assert.That(resolved, Is.TypeOf<FromFileDeconvolutionParameters>());
            Assert.That(((FromFileDeconvolutionParameters)resolved).MinAssumedChargeState, Is.EqualTo(3));
            Assert.That(((FromFileDeconvolutionParameters)resolved).MaxAssumedChargeState, Is.EqualTo(20));

            // Materialization via SetAllFileSpecificCommonParams (the CLI path) also works
            var fileSpecific = MetaMorpheusTask.SetAllFileSpecificCommonParams(
                loadedTask.CommonParameters, null, rawPath);
            Assert.That(fileSpecific.PrecursorDeconvolutionParameters, Is.TypeOf<FromFileDeconvolutionParameters>());
            Assert.That(((FromFileDeconvolutionParameters)fileSpecific.PrecursorDeconvolutionParameters).MinAssumedChargeState, Is.EqualTo(3));

            // SourceMapPath exclusion from equality is verified by SearchFeatureFileMap_Equality_IgnoresSourceMapPath
            // This test proves the functional contract: SourceMapPath is audit-only, embedded entries drive materialization
        }
        finally
        {
            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir, true);
        }
    }

    #endregion
}
