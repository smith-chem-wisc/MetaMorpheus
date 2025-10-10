using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Threading.Tasks;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using NUnit.Framework;
using Omics;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;
using TaskLayer;
using Transcriptomics;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.Transcriptomics;

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
        var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", databaseFileName);
        var dbForTask = new List<DbForTask> { new DbForTask(dbPath, false) };
        var commonParameters = new CommonParameters(digestionParams: new RnaDigestionParams());

        var loader = new DatabaseLoadingEngine(commonParameters, [], [], dbForTask, "TestTaskId", DecoyType.None);
        var results = (DatabaseLoadingEngineResults)loader.Run()!;
        var bioPolymers = results.BioPolymers;

        Assert.That(dbForTask[0].BioPolymerCount, Is.EqualTo(1));   
        Assert.That(dbForTask[0].TargetCount, Is.EqualTo(1));
        Assert.That(dbForTask[0].DecoyCount, Is.EqualTo(0));

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
        var commonParameters = new CommonParameters(digestionParams: new RnaDigestionParams());
        string baseSequence = "GUACUGCCUCUAGUGAAGCA";

        // Create two instances of the same database, one as a contaminant, one as a target. 
        var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "20mer1.fasta");
        var dbsForTask = new List<DbForTask> { new DbForTask(dbPath, false), new DbForTask(dbPath, true) };
        var loader = new DatabaseLoadingEngine(commonParameters, [], [], dbsForTask, "TestTaskId", DecoyType.None, tcAmbiguity: type);
        var results = (DatabaseLoadingEngineResults)loader.Run()!;
        var bioPolymers = results.BioPolymers;

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

        var loader = new DatabaseLoadingEngine(new(), [], [], [], "TestTaskId", DecoyType.None, tcAmbiguity: TargetContaminantAmbiguity.RenameProtein);
        var exception = Assert.Throws<ArgumentException>(() =>
        {
            loader.SanitizeBioPolymerDatabase(notImplementedBioPolymers, TargetContaminantAmbiguity.RenameProtein);
        });

        Assert.That(exception, Is.Not.Null);
        Assert.That(exception.Message, Does.Contain("Database sanitization assumed BioPolymer was a protein when it was Test.Transcriptomics.TestBioPolymer"));
    }

    [Test]
    public void TwoTruncationsAndSequenceVariant_DbLoading()
    {
        GlobalVariables.AnalyteType = AnalyteType.Oligo;
        var task = new SearchTask();
        var commonParameters = new CommonParameters(digestionParams: new RnaDigestionParams("RNase T1"));
        string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "TruncationAndVariantMods.xml");
        var dbsForTask = new List<DbForTask> { new DbForTask(dbPath, false) };

        // Use reflection to access the protected LoadModifications method
        var loadModificationsMethod = typeof(MetaMorpheusTask).GetMethod("LoadModifications", BindingFlags.NonPublic | BindingFlags.Instance);
        object[] modArgs = new object[]
        {
            "TestTaskId",
            null, // variableModifications (out)
            null, // fixedModifications (out)
            null  // localizableModificationTypes (out)
        };
        loadModificationsMethod.Invoke(task, modArgs);
        var localizableModificationTypes = (List<string>)modArgs[3];

        var loader = new DatabaseLoadingEngine(commonParameters, [], [], dbsForTask, "TestTaskId", DecoyType.Reverse, true, localizableModificationTypes);
        var results = (DatabaseLoadingEngineResults)loader.Run()!;
        var rna = results.BioPolymers;

        Assert.That(rna.All(p => p.SequenceVariations.Count == 1));
        Assert.That(rna.All(p => p.OriginalNonVariantModifications.Count == 2));

        List<IBioPolymer> targets = rna.Where(p => p.IsDecoy == false).ToList();
        IBioPolymer variantTarget = targets.First(p => p.AppliedSequenceVariations.Count >= 1);
        IBioPolymer nonVariantTarget = targets.First(p => p.AppliedSequenceVariations.Count == 0);

        Assert.That(variantTarget.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
        Assert.That(nonVariantTarget.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));

        List<IBioPolymer> decoys = rna.Where(p => p.IsDecoy).ToList();
        IBioPolymer variantDecoy = decoys.First(p => p.AppliedSequenceVariations.Count >= 1);
        IBioPolymer nonVariantDecoy = decoys.First(p => p.AppliedSequenceVariations.Count == 0);

        Assert.That(variantDecoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
        Assert.That(nonVariantDecoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
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