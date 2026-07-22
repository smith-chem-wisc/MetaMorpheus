using System;
using System.IO;
using NUnit.Framework;
using EngineLayer.Deconvolution.FeatureFileMapping;

namespace Test
{
    [TestFixture]
    public class FeatureFileMapStoreTests
    {
        private string tempDir;
        private string testStorePath;

        [SetUp]
        public void Setup()
        {
            tempDir = Path.Combine(Path.GetTempPath(), Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(tempDir);
            testStorePath = Path.Combine(tempDir, "feature-maps.toml");
        }

        [TearDown]
        public void Teardown()
        {
            if (Directory.Exists(tempDir))
            {
                Directory.Delete(tempDir, true);
            }
        }

        [Test]
        public void StoreService_FirstRun_ReturnsEmptyMap()
        {
            // Assert file does not exist
            Assert.That(File.Exists(testStorePath), Is.False);

            // Act
            var map = FeatureFileMapStore.Load(testStorePath);

            // Assert
            Assert.That(map, Is.Not.Null);
            Assert.That(map.Count, Is.EqualTo(0));
            Assert.That(map.Conditions.Count, Is.EqualTo(0));
            Assert.That(map.FilePath, Is.EqualTo(Path.GetFullPath(testStorePath)));
        }

        [Test]
        public void StoreService_SaveAndReload_PreservesData()
        {
            var map = new FeatureFileMap();
            var condition = new FeatureFileMapCondition { ConditionKey = "Cond1", DisplayName = "Condition 1" };
            map.AddOrReplaceCondition(condition);

            var entry = new FeatureSpectraEntry { MassSpecFilePath = "raw1.raw", MassSpecFileName = "raw1.raw" };
            entry.FeatureFiles.Add(new FeatureFileConditionEntry { ConditionKey = "Cond1", FeatureFilePath = "feature1.tsv", FeatureFileName = "feature1.tsv" });
            map.AddOrReplaceSpectraEntry(entry);

            // Act
            FeatureFileMapStore.Save(map, testStorePath);
            Assert.That(File.Exists(testStorePath), Is.True);

            var reloaded = FeatureFileMapStore.Load(testStorePath);

            // Assert
            Assert.That(reloaded.Conditions.Count, Is.EqualTo(1));
            Assert.That(reloaded.Conditions[0].ConditionKey, Is.EqualTo("Cond1"));

            Assert.That(reloaded.Count, Is.EqualTo(1));
            Assert.That(reloaded.TryGetSpectraEntry(Path.GetFullPath("raw1.raw"), out var reloadedEntry), Is.True);
            Assert.That(reloadedEntry.FeatureFiles.Count, Is.EqualTo(1));
            Assert.That(reloadedEntry.FeatureFiles[0].FeatureFilePath, Is.EqualTo(Path.GetFullPath("feature1.tsv")));
        }

        [Test]
        public void StoreService_InvalidToml_ThrowsFeatureMappingException()
        {
            File.WriteAllText(testStorePath, "This is not valid TOML = { [ ]");

            Assert.Throws<FeatureMappingException>(() => FeatureFileMapStore.Load(testStorePath));
        }

        [Test]
        public void StoreService_PathNormalization_CollapsesCasingDifferences()
        {
            var map = new FeatureFileMap();
            
            // Note: Windows path casing normalization
            string upperPath = Path.Combine(tempDir, "RAW1.RAW");
            string lowerPath = Path.Combine(tempDir, "raw1.raw");

            var entry = new FeatureSpectraEntry { MassSpecFilePath = upperPath, MassSpecFileName = "RAW1.RAW" };
            map.AddOrReplaceSpectraEntry(entry);

            FeatureFileMapStore.Save(map, testStorePath);

            var reloaded = FeatureFileMapStore.Load(testStorePath);

            // Act & Assert
            // TryGetSpectraEntry uses MatchesRawFile which internally calls NormalizePath
            Assert.That(reloaded.TryGetSpectraEntry(lowerPath, out var retrievedEntry), Is.True);
            Assert.That(Path.GetFullPath(retrievedEntry.MassSpecFilePath), Is.EqualTo(Path.GetFullPath(upperPath)));
        }

        [Test]
        public void StoreService_AtomicRewrite_HandlesExistingFile()
        {
            var map = new FeatureFileMap();
            map.AddOrReplaceCondition(new FeatureFileMapCondition { ConditionKey = "Cond1" });
            FeatureFileMapStore.Save(map, testStorePath);
            
            Assert.That(File.Exists(testStorePath), Is.True);

            // Write again
            map.AddOrReplaceCondition(new FeatureFileMapCondition { ConditionKey = "Cond2" });
            FeatureFileMapStore.Save(map, testStorePath);

            var reloaded = FeatureFileMapStore.Load(testStorePath);
            Assert.That(reloaded.Conditions.Count, Is.EqualTo(2));
        }
    }
}
