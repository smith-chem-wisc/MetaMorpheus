#nullable enable
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer;
using EngineLayer.DiaLibrarySearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

namespace Test.DiaLibrarySearch;

/// <summary>
/// Protein quantification for the DIA search goes through mzLib's QuantificationEngine, the same engine MetaMorpheus
/// uses for isobaric quantification. Each identified precursor becomes one spectral match whose single intensity is
/// its quantity, attributed to exactly one peptide: the one parsimony assigned to its protein group.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaProteinQuantificationTests
{
    // A has EEEEER of its own, C has VVVVVK and TTTTTR; LLLLLK is in A and B. Parsimony drops B, so among the groups
    // reported LLLLLK is A's alone: group-level uniqueness (mzLib #1463, MetaMorpheus 1.1.13, DEF-PROT-UNIQUE v3.14).
    private static readonly List<Protein> Targets =
    [
        new("LLLLLKEEEEER", "A"),
        new("LLLLLKGGGGGR", "B"),
        new("VVVVVKTTTTTR", "C"),
    ];

    private static readonly CommonParameters Parameters = new(digestionParams: new DigestionParams(minPeptideLength: 5));

    private static List<Protein> Database() =>
        Targets.Concat(DecoyProteinGenerator.GenerateDecoys(Targets, DecoyType.Reverse)).ToList();

    private static DiaPrecursorMatch Precursor(int index, string sequence, int charge, double quantity, double qValue = 0.001) =>
        new(index, sequence, charge, 500, false, new Irt(0), new RtMinutes(10), new Irt(0), 5, [5.0], qValue, quantity);

    private static readonly List<DiaPrecursorMatch> Run1 =
    [
        Precursor(0, "EEEEER", 2, 100),
        Precursor(1, "EEEEER", 3, 50),
        Precursor(2, "LLLLLK", 2, 1000),
        Precursor(3, "VVVVVK", 2, 30),
        Precursor(4, "TTTTTR", 2, 70),
    ];

    private static readonly List<DiaPrecursorMatch> Run2 =
    [
        Precursor(0, "EEEEER", 2, 300),
        Precursor(3, "VVVVVK", 2, 10),
    ];

    private static DiaProteinInferenceResult Infer(string file, IReadOnlyList<DiaPrecursorMatch> precursors) =>
        DiaProteinInference.Infer(precursors.GroupBy(p => p.FullSequence)
            .Select(g => new DiaPeptideMatch(g.Key, false, g.First().PrecursorIndex, 5, g.Min(p => p.QValue))).ToList(),
            Database(), Parameters, file, [], []);

    [Test]
    public void APrecursorIsOneMatchWithItsQuantityAsItsOnlyIntensity()
    {
        var inference = Infer(@"C:\data\run1.raw", Run1);

        var quantified = DiaProteinQuantification.QuantifiedPrecursors(@"C:\data\run1.raw", Run1, inference, 0.01);

        var e2 = quantified.Single(q => q.FullSequence == "EEEEER" && q.OneBasedScanNumber == 1);
        Assert.That(e2.Intensities, Is.EqualTo(new[] { 100.0 }));
        Assert.That(e2.FullFilePath, Is.EqualTo(@"C:\data\run1.raw"));
        Assert.That(e2.Accession, Is.EqualTo("A"));
        Assert.That(e2.IsDecoy, Is.False);
        Assert.That(e2.GetIdentifiedBioPolymersWithSetMods().Single(),
            Is.SameAs(inference.ProteinGroups.Single(g => g.ProteinGroupName == "A").UniquePeptides.Single(p => p.FullSequence == "EEEEER")));
    }

    /// <summary>A shared peptide is still attributed to one peptide, never left ambiguous for the engine to drop.</summary>
    [Test]
    public void ASharedPeptideIsAttributedToOnePeptide()
    {
        var inference = Infer(@"C:\data\run1.raw", Run1);

        var shared = DiaProteinQuantification.QuantifiedPrecursors(@"C:\data\run1.raw", Run1, inference, 0.01)
            .Single(q => q.FullSequence == "LLLLLK");

        Assert.That(shared.GetIdentifiedBioPolymersWithSetMods().Count(), Is.EqualTo(1));
    }

    /// <summary>
    /// Only reported precursors are quantified: a target passing the threshold, whose peptide inference mapped, with a
    /// measured quantity. An unmeasured quantity is null, never NaN or zero, so the engine leaves it out.
    /// </summary>
    [Test]
    public void OnlyReportedAndMeasuredPrecursorsAreQuantified()
    {
        var precursors = new List<DiaPrecursorMatch>
        {
            Precursor(0, "EEEEER", 2, 100),
            Precursor(1, "VVVVVK", 2, 30, qValue: 0.5),
            Precursor(2, "TTTTTR", 2, double.NaN),
            Precursor(3, "WWWWWK", 2, 40),
        };
        var inference = Infer(@"C:\data\run1.raw", precursors);

        var quantified = DiaProteinQuantification.QuantifiedPrecursors(@"C:\data\run1.raw", precursors, inference, 0.01);

        Assert.That(quantified.Select(q => q.FullSequence), Is.EquivalentTo(new[] { "EEEEER", "TTTTTR" }));
        Assert.That(quantified.Single(q => q.FullSequence == "TTTTTR").Intensities, Is.Null);
    }

    /// <summary>
    /// Precursors sum to peptides and unique peptides sum to protein groups, per run, with nothing normalised: the
    /// conservative strategies MetaMorpheus already uses for isobaric quantification.
    /// </summary>
    [Test]
    public void ProteinsAreTheSumOfTheirUniquePeptidesInEachRun()
    {
        var inference = Infer(@"C:\data\run1.raw", Run1);
        var quantified = DiaProteinQuantification.QuantifiedPrecursors(@"C:\data\run1.raw", Run1, inference, 0.01)
            .Concat(DiaProteinQuantification.QuantifiedPrecursors(@"C:\data\run2.raw", Run2, inference, 0.01))
            .ToList();

        var results = DiaProteinQuantification.Quantify(quantified, inference.ProteinGroups, outputDirectory: null);

        Assert.That(results.Success, Is.True, results.Summary);
        var a = results.ProteinIntensities.Single(p => p.Key.BioPolymerGroupName == "A").Value;
        var c = results.ProteinIntensities.Single(p => p.Key.BioPolymerGroupName == "C").Value;
        // Run 1: EEEEER 100 + 50, and LLLLLK 1000, unique to A once B is dropped. Run 2: EEEEER 300.
        Assert.That(a.OrderBy(kv => ((SpectraFileInfo)kv.Key).FullFilePathWithExtension).Select(kv => kv.Value), Is.EqualTo(new[] { 1150.0, 300.0 }));
        Assert.That(c.OrderBy(kv => ((SpectraFileInfo)kv.Key).FullFilePathWithExtension).Select(kv => kv.Value), Is.EqualTo(new[] { 100.0, 10.0 }));
        Assert.That(results.AmbiguousSpectralMatchesExcluded, Is.EqualTo(0));
        Assert.That(results.WrittenFiles, Is.Empty);
    }
}
