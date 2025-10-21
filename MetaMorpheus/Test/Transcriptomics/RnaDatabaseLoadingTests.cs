using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Threading.Tasks;
using EngineLayer;
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
        var bioPolymers = new SearchTask().LoadBioPolymers(
            "TestTaskId",
            dbForTask,
            true, 
            DecoyType.None,
            [], 
            commonParameters
        );

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
        var bioPolymers = new SearchTask().LoadBioPolymers(
            "TestTaskId",
            dbsForTask,
            true,
            DecoyType.None,
            [],
            commonParameters
        );

        // Use reflection to access protected SanitizeBiopolymers method
        var sanitizeBioPolymersMethod = typeof(MetaMorpheusTask).GetMethod("SanitizeBioPolymerDatabase",
            BindingFlags.NonPublic | BindingFlags.Static);
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
    [Test]
    public void TwoTruncationsAndSequenceVariant_DbLoading()
    {
        // Purpose of this test:
        // - Verify RNA DB loading from an XML that encodes two truncations and one sequence variant.
        // - Confirm that "nucleotide substitution" PTMs from XML are converted into true SequenceVariations.
        // - Ensure two truncation-related non-variant modifications are present on each entry.
        // - Confirm the loader emits both canonical and variant-applied isoforms (for targets and decoys).
        // - Validate that applying the variant reduces the number of valid localized mods for that isoform (2 -> 1),
        //   because per-isoform mods are revalidated after the sequence change.
        //
        // Relevant mechanics from ProteinDbLoader:
        // - LoadProteinXML reads Protein XML entries and builds biopolymers.
        // - If OneBasedPossibleLocalizedModifications contains a "nucleotide substitution" type, the loader calls
        //   ConvertNucleotideSubstitutionModificationsToSequenceVariants<T>() to promote it to SequenceVariations.
        // - Variants are expanded via GetVariantBioPolymers(...), yielding canonical (no AppliedSequenceVariations) and
        //   variant-applied isoforms (AppliedSequenceVariations.Count >= 1).
        // - Per-isoform OneBasedPossibleLocalizedModifications are revalidated on the new sequence; incompatible mods are dropped.
        // - Decoys are generated (here: DecoyType.Reverse) after target loading and mirror the variant/truncation logic.

        // Step 1: Configure RNA analyte handling so the pipeline treats entries as nucleic-acid biopolymers (RNA).
        GlobalVariables.AnalyteType = AnalyteType.Oligo;

        // Step 2: Prepare a SearchTask and CommonParameters with RNase T1 digestion (not directly exercised here,
        // but consistent with normal configuration).
        // Request at least one variant-applied isoform to be emitted
        var commonParameters = new CommonParameters(
            digestionParams: new RnaDigestionParams("RNase T1"),
            maxSequenceVariantsPerIsoform: 1,
            minAlleleDepth: 0,
            maxSequenceVariantIsoforms: 2
        );

        var task = new SearchTask
        {
            CommonParameters = commonParameters
        };

        // Step 3: Use TruncationAndVariantMods.xml which encodes exactly one sequence variant and two truncations.
        string dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "TruncationAndVariantMods.xml");
        var dbsForTask = new List<DbForTask> { new DbForTask(dbPath, false) };

        // Step 4: Load localizable modification types the same way the pipeline would (via reflection on MetaMorpheusTask.LoadModifications).
        // Passing these into LoadBioPolymers ensures consistent classification/filtering of modifications during DB load.
        var loadModificationsMethod = typeof(MetaMorpheusTask).GetMethod("LoadModifications", BindingFlags.NonPublic | BindingFlags.Instance);
        object[] modArgs = new object[] { "TestTaskId", null, null, null };
        loadModificationsMethod.Invoke(task, modArgs);
        var localizableModificationTypes = (List<string>)modArgs[3];

        // Additional NUnit assertions (sanity): localizable modification types should be available.
        Assert.Multiple(() =>
        {
            Assert.That(localizableModificationTypes, Is.Not.Null, "Localizable modification types should not be null.");
            Assert.That(localizableModificationTypes, Is.Not.Empty, "Expected at least one localizable modification type.");
        });

        // Step 5: Load RNA biopolymers with reverse decoys. Internally, XML "nucleotide substitution" PTMs are
        // converted to SequenceVariations and variants are expanded into canonical + variant-applied isoforms.
        var rna = task.LoadBioPolymers(
            "TestTaskId",
            dbsForTask,
            true,
            DecoyType.Reverse,
            localizableModificationTypes,
            commonParameters
        );

        // Step 6: High-level expectations across all returned biopolymers (targets + decoys).
        // - Exactly one potential SequenceVariation per entry.
        // - Exactly two original (non-variant) modifications per entry (the two truncations).
        Assert.Multiple(() =>
        {
            Assert.That(rna, Is.Not.Null.And.Not.Empty, "DB load should produce at least one biopolymer.");
            Assert.That(rna.All(p => p.SequenceVariations.Count == 1), Is.True,
                "Each entry should carry exactly one possible SequenceVariation.");
            Assert.That(rna.All(p => p.OriginalNonVariantModifications.Count == 2), Is.True,
                "Each entry should carry exactly two original (non-variant) modifications (two truncations).");
            Assert.That(rna.Any(p => p.IsDecoy), Is.True, "Decoys should be present when DecoyType.Reverse is used.");
            Assert.That(rna.Any(p => !p.IsDecoy), Is.True, "Targets should be present when searchTarget is true.");
        });

        // Step 7: Partition targets and verify we have both canonical and variant-applied isoforms.
        List<IBioPolymer> targets = rna.Where(p => p.IsDecoy == false).ToList();

        // Assert first to get clear failure messages before First(...) could throw.
        Assert.Multiple(() =>
        {
            Assert.That(targets.Any(p => p.AppliedSequenceVariations.Count >= 1), Is.True,
                "Expected at least one variant-applied target isoform.");
            Assert.That(targets.Any(p => p.AppliedSequenceVariations.Count == 0), Is.True,
                "Expected at least one canonical (non-variant) target isoform.");
        });

        IBioPolymer variantTarget = targets.First(p => p.AppliedSequenceVariations.Count >= 1);
        IBioPolymer nonVariantTarget = targets.First(p => p.AppliedSequenceVariations.Count == 0);

        // After variant application, per-isoform OneBasedPossibleLocalizedModifications are revalidated against the altered sequence,
        // so the variant isoform should drop one incompatible localized mod (2 -> 1).
        Assert.That(variantTarget.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
        Assert.That(nonVariantTarget.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));

        // Additional supportive asserts for targets
        Assert.Multiple(() =>
        {
            // Canonical target should have no applied variants; variant target should have at least one applied variant.
            Assert.That(nonVariantTarget.AppliedSequenceVariations, Is.Empty, "Canonical target should have no applied variants.");
            Assert.That(variantTarget.AppliedSequenceVariations.Count, Is.GreaterThanOrEqualTo(1),
                "Variant target should have one or more applied variants.");

            // The base sequence should change between canonical and variant-applied isoforms.
            Assert.That(variantTarget.BaseSequence, Is.Not.EqualTo(nonVariantTarget.BaseSequence),
                "Variant target base sequence should differ from the canonical base sequence.");

            // The non-variant mod set should be consistent with earlier global assertion (2).
            Assert.That(nonVariantTarget.OriginalNonVariantModifications.Count, Is.EqualTo(2));

            // Validate that variant isoform's localized mod sites are a subset of the canonical's (due to dropping incompatibles).
            var variantKeys = variantTarget.OneBasedPossibleLocalizedModifications.Keys.ToHashSet();
            var canonicalKeys = nonVariantTarget.OneBasedPossibleLocalizedModifications.Keys.ToHashSet();
            Assert.That(variantKeys.IsSubsetOf(canonicalKeys), Is.True,
                "Variant isoform localized modification sites should be a subset of canonical sites.");
        });

        // Step 8: Partition decoys and verify the same canonical/variant-applied behavior is mirrored.
        List<IBioPolymer> decoys = rna.Where(p => p.IsDecoy).ToList();

        // Assert first to get clear failure messages before First(...) could throw.
        Assert.Multiple(() =>
        {
            Assert.That(decoys.Any(p => p.AppliedSequenceVariations.Count >= 1), Is.True,
                "Expected at least one variant-applied decoy isoform.");
            Assert.That(decoys.Any(p => p.AppliedSequenceVariations.Count == 0), Is.True,
                "Expected at least one canonical (non-variant) decoy isoform.");
        });

        IBioPolymer variantDecoy = decoys.First(p => p.AppliedSequenceVariations.Count >= 1);
        IBioPolymer nonVariantDecoy = decoys.First(p => p.AppliedSequenceVariations.Count == 0);

        // Decoys should preserve the same logical relationship between variants and truncation mods after reversal.
        Assert.That(variantDecoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
        Assert.That(nonVariantDecoy.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));

        // Additional supportive asserts for decoys
        Assert.Multiple(() =>
        {
            // Canonical decoy vs variant decoy, applied variant presence.
            Assert.That(nonVariantDecoy.AppliedSequenceVariations, Is.Empty, "Canonical decoy should have no applied variants.");
            Assert.That(variantDecoy.AppliedSequenceVariations.Count, Is.GreaterThanOrEqualTo(1),
                "Variant decoy should have one or more applied variants.");

            // Decoy variant vs canonical base sequences should differ.
            Assert.That(variantDecoy.BaseSequence, Is.Not.EqualTo(nonVariantDecoy.BaseSequence),
                "Variant decoy base sequence should differ from the canonical decoy base sequence.");

            // Validate that variant decoy's localized mod sites are a subset of the canonical decoy's.
            var variantDecoyKeys = variantDecoy.OneBasedPossibleLocalizedModifications.Keys.ToHashSet();
            var canonicalDecoyKeys = nonVariantDecoy.OneBasedPossibleLocalizedModifications.Keys.ToHashSet();
            Assert.That(variantDecoyKeys.IsSubsetOf(canonicalDecoyKeys), Is.True,
                "Variant decoy localized modification sites should be a subset of canonical decoy sites.");
        });
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