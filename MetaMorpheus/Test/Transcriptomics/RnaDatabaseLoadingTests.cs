using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Drawing.Text;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Reflection.Metadata;
using EngineLayer;
using NUnit.Framework;
using Omics;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;
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

    [Test]
    public static void DatabaseSanitization_ThrowsWhenNewBioPolymerIntroduced()
    {
        List<IBioPolymer> notImplementedBioPolymers = new()
        {
            new TestBioPolymer() {Accession= "Testing"},
            new TestBioPolymer() {Accession= "Testing"},
            new TestBioPolymer() {Accession= "TestingA"},
        };

        // Use reflection to access protected SanitizeBiopolymers method
        var sanitizeBioPolymersMethod = typeof(MetaMorpheusTask).GetMethod("SanitizeBioPolymerDatabase",
            BindingFlags.NonPublic | BindingFlags.Static);
        var genericMethod = sanitizeBioPolymersMethod.MakeGenericMethod(typeof(IBioPolymer));

        var exception = Assert.Throws<TargetInvocationException>(() =>
        {
            genericMethod.Invoke(null, new object[] { notImplementedBioPolymers, TargetContaminantAmbiguity.RenameProtein });
        });

        Assert.That(exception, Is.Not.Null);
        Assert.That(exception.InnerException, Is.Not.Null);
        Assert.That(exception.InnerException, Is.TypeOf<ArgumentException>());
        Assert.That(exception.InnerException.Message, Does.Contain("Database sanitization assumed BioPolymer was a protein when it was Test.Transcriptomics.TestBioPolymer"));
    }
}


[ExcludeFromCodeCoverage]
class TestBioPolymer : IBioPolymer
{
    private string _baseSequence;
    private IDictionary<int, List<Modification>> _oneBasedPossibleLocalizedModifications;
    private string _baseSequence1;
    private IDictionary<int, List<Modification>> _oneBasedPossibleLocalizedModifications1;

    public TBioPolymerType CreateVariant<TBioPolymerType>(string variantBaseSequence, TBioPolymerType original,
        IEnumerable<SequenceVariation> appliedSequenceVariants, IEnumerable<TruncationProduct> applicableProteolysisProducts, IDictionary<int, List<Modification>> oneBasedModifications,
        string sampleNameForVariants) where TBioPolymerType : IHasSequenceVariants
    {
        throw new NotImplementedException();
    }

    /// <summary>
    /// Basic Constructor that creates a TestBioPolymer object with a default sequence and no modifications.
    /// We can fill in details of this class as new testing situations arise. 
    /// </summary>
    /// <param name="baseSequence"></param>
    /// <param name="baseSequence1"></param>
    /// <param name="oneBasedPossibleLocalizedModifications"></param>
    /// <param name="oneBasedPossibleLocalizedModifications1"></param>
    public TestBioPolymer(string baseSequence = "GUACUGCCUCUAGUGAAGCA", string baseSequence1 = "GUACUGCCUCUAGUGAAGCA",
        IDictionary<int, List<Modification>> oneBasedPossibleLocalizedModifications = null,
        IDictionary<int, List<Modification>> oneBasedPossibleLocalizedModifications1 = null)
    {
        _baseSequence = baseSequence;
        _baseSequence1 = baseSequence1;
        _oneBasedPossibleLocalizedModifications = oneBasedPossibleLocalizedModifications ?? new Dictionary<int, List<Modification>>();
        _oneBasedPossibleLocalizedModifications1 = oneBasedPossibleLocalizedModifications1 ?? new Dictionary<int, List<Modification>>();
        Length = baseSequence.Length;
        DatabaseFilePath = string.Empty;
        IsDecoy = false;
        IsContaminant = false;
        Organism = string.Empty;
        Accession = string.Empty;
        GeneNames = new List<Tuple<string, string>>();
        Name = string.Empty;
        FullName = string.Empty;
        SampleNameForVariants = string.Empty;
        OriginalNonVariantModifications = new Dictionary<int, List<Modification>>();
        ConsensusVariant = null;
        AppliedSequenceVariations = new List<SequenceVariation>();
        SequenceVariations = new List<SequenceVariation>();
        TruncationProducts = new List<TruncationProduct>();
    }

    string IBioPolymer.BaseSequence => _baseSequence1;

    public int Length { get; }
    public string DatabaseFilePath { get; }
    public bool IsDecoy { get; }
    public bool IsContaminant { get; }
    public string Organism { get; }
    public string Accession { get; set; }
    public List<Tuple<string, string>> GeneNames { get; }

    IDictionary<int, List<Modification>> IBioPolymer.OneBasedPossibleLocalizedModifications => _oneBasedPossibleLocalizedModifications1;

    public IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams digestionParams, List<Modification> allKnownFixedModifications, List<Modification> variableModifications,
        List<SilacLabel> silacLabels = null, (SilacLabel startLabel, SilacLabel endLabel)? turnoverLabels = null,
        bool topDownTruncationSearch = false)
    {
        throw new NotImplementedException();
    }

    public IBioPolymer CloneWithNewSequenceAndMods(string newBaseSequence, IDictionary<int, List<Modification>> newMods)
    {
        throw new NotImplementedException();
    }

    public string Name { get; }
    public string FullName { get; }

    string IHasSequenceVariants.BaseSequence => _baseSequence;

    public string SampleNameForVariants { get; }

    IDictionary<int, List<Modification>> IHasSequenceVariants.OneBasedPossibleLocalizedModifications => _oneBasedPossibleLocalizedModifications;

    public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; }
    public IBioPolymer ConsensusVariant { get; }
    public List<SequenceVariation> AppliedSequenceVariations { get; }
    public List<SequenceVariation> SequenceVariations { get; }
    public List<TruncationProduct> TruncationProducts { get; }
}