using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection.Metadata;
using EngineLayer;
using NUnit.Framework;
using Omics;
using TaskLayer;
using Transcriptomics;
using UsefulProteomicsDatabases;

namespace Test.Transcriptomics;

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

    [Test]
    [TestCase(TargetContaminantAmbiguity.RenameProtein, new[] { "20mer1_D1", "20mer1_D2" }, new[] {false, true})]
    [TestCase(TargetContaminantAmbiguity.RemoveTarget, new[] { "20mer1" }, new[] { true })]
    [TestCase(TargetContaminantAmbiguity.RemoveContaminant, new[] { "20mer1" }, new[] { false })]
    public static void DbReadingHandleAccessionCollisions(TargetContaminantAmbiguity type, string[] expectedAccessions, bool[] expectedIsContaminant)
    {
        GlobalVariables.AnalyteType = AnalyteType.Oligo;
        var task = new SearchTask();
        var commonParameters = new CommonParameters();
        string baseSequence = "GUACUGCCUCUAGUGAAGCA";

        // Create two instances of the same database, one as a contaminant, one as a target. 
        var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "20mer1.fasta");
        var dbsForTask = new List<DbForTask> { new DbForTask(dbPath, false), new DbForTask(dbPath, true) };

        // Use reflection to access the protected LoadBioPolymers method
        var loadBioPolymersMethod = typeof(MetaMorpheusTask).GetMethod("LoadBioPolymers", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        if (loadBioPolymersMethod == null)
        {
            Assert.Fail("Failed to find the LoadBioPolymers method via reflection.");
        }

        var bioPolymers = (List<IBioPolymer>)loadBioPolymersMethod.Invoke(task, new object[]
        {
            "TestTaskId",
            dbsForTask,
            true, // searchTarget
            DecoyType.None,
            new List<string>(), // localizeableModificationTypes
            commonParameters
        });

        // Use reflection to access protected SanitizeBiopolymers method
        var sanitizeBioPolymersMethod = typeof(MetaMorpheusTask).GetMethod("SanitizeBioPolymerDatabase",
            System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);
        var genericMethod = sanitizeBioPolymersMethod.MakeGenericMethod(typeof(IBioPolymer));
        genericMethod.Invoke(null, new object[] { bioPolymers, type });

        Assert.That(bioPolymers.Count, Is.EqualTo(expectedAccessions.Length));
        Assert.That(bioPolymers[0].Accession, Is.EqualTo(expectedAccessions[0]));
        Assert.That(bioPolymers[0].IsContaminant, Is.EqualTo(expectedIsContaminant[0]));

        for (int i = 0; i < bioPolymers.Count; i++)
        {
            Assert.That(bioPolymers[i].BaseSequence, Is.EqualTo(baseSequence));
            Assert.That(bioPolymers[i].Accession, Is.EqualTo(expectedAccessions[i]));
            Assert.That(bioPolymers[i].IsContaminant, Is.EqualTo(expectedIsContaminant[i]));
        }
    }
}