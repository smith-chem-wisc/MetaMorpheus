using NUnit.Framework;
using EngineLayer;
using Omics.Modifications;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Transcriptomics;
using UsefulProteomicsDatabases.Transcriptomics;
using UsefulProteomicsDatabases;
using Omics;
using TaskLayer;

namespace Test.Transcriptomics;

[TestFixture]
public class RnaModLoadingTests
{
    [Test]
    public static void ModsLoadWithoutCrash()
    {
        // This test is to ensure that the mods load without crashing.
        // It does not check for correctness of the mods.

        Assert.That(GlobalVariables.AllRnaModsKnown, Is.Not.Null);
        Assert.That(GlobalVariables.AllRnaModsKnown, Is.Not.Empty);

        Assert.That(GlobalVariables.AllRnaModTypesKnown, Is.Not.Null);
        Assert.That(GlobalVariables.AllRnaModTypesKnown, Is.Not.Empty);

        Assert.That(GlobalVariables.AllRnaModsKnownDictionary, Is.Not.Null);
        Assert.That(GlobalVariables.AllRnaModsKnownDictionary, Is.Not.Empty);
    }
}

[TestFixture]
public class RnaDatabaseLoadingTests
{
    [Test]
    [TestCase("20mer1.fasta")]
    [TestCase("20mer1.fasta.gz")]
    [TestCase("20mer1.xml")]
    [TestCase("20mer1.xml.gz")]
    public static void TestDbReadingDifferentExtensions(string databaseFileName)
    {
        GlobalVariables.AnalyteType = AnalyteType.Oligo;
        List<IBioPolymer> bioPolymers = new List<IBioPolymer>();
        var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", databaseFileName);

        // Use reflection to access the protected LoadBioPolymers method
        var loadBioPolymersMethod = typeof(MetaMorpheusTask).GetMethod("LoadBioPolymers", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        if (loadBioPolymersMethod == null)
        {
            Assert.Fail("Failed to find the LoadBioPolymers method via reflection.");
        }

        var task = new SearchTask();
        var dbForTask = new List<DbForTask> { new DbForTask(dbPath, false) };
        var commonParameters = new CommonParameters();

        // Invoke the method
        bioPolymers = (List<IBioPolymer>)loadBioPolymersMethod.Invoke(task, new object[]
        {
            "TestTaskId",
            dbForTask,
            true, // searchTarget
            DecoyType.None,
            new List<string>(), // localizeableModificationTypes
            commonParameters
        });

        Assert.That(bioPolymers![0], Is.TypeOf<RNA>());
        Assert.That(bioPolymers.Count, Is.EqualTo(1));
        Assert.That(bioPolymers.First().BaseSequence, Is.EqualTo("GUACUGCCUCUAGUGAAGCA"));
    }
}
