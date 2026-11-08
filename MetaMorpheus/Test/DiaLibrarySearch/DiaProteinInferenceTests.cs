#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer;
using EngineLayer.DiaLibrarySearch;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

namespace Test.DiaLibrarySearch;

/// <summary>
/// Protein inference for the DIA search reuses MetaMorpheus's DDA engines unchanged (ProteinParsimonyEngine, then
/// ProteinScoringAndFdrEngine). The adapter's only job is to present each DIA peptide as the spectral match those
/// engines expect: mapped to every protein that yields it, carrying its peptide q-value and a positive score.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaProteinInferenceTests
{
    // Three proteins. LLLLLK is shared by A and B; EEEEER is unique to A, GGGGGR to B, and VVVVVK and TTTTTR to C.
    private static readonly List<Protein> Targets =
    [
        new("LLLLLKEEEEER", "A"),
        new("LLLLLKGGGGGR", "B"),
        new("VVVVVKTTTTTR", "C"),
    ];

    private static readonly CommonParameters Parameters = new(digestionParams: new DigestionParams(minPeptideLength: 5));

    private static List<Protein> TargetsAndDecoys() =>
        Targets.Concat(DecoyProteinGenerator.GenerateDecoys(Targets, DecoyType.Reverse)).ToList();

    private static Dictionary<string, List<PeptideWithSetModifications>> Index() =>
        DiaProteinInference.IndexPeptides(TargetsAndDecoys(), Parameters, [], []);

    private static DiaPeptideMatch Peptide(string sequence, double qValue, double score = 5, bool isDecoy = false, int index = 0) =>
        new(sequence, isDecoy, index, score, qValue);

    private static string DecoySequenceOf(string accession) =>
        Index().Where(e => e.Value.All(p => p.Parent.Accession == "DECOY_" + accession)).Select(e => e.Key)
            .OrderBy(s => s, StringComparer.Ordinal).First();

    [Test]
    public void APeptideMapsToEveryProteinThatYieldsIt()
    {
        var index = Index();

        Assert.That(index["LLLLLK"].Select(p => p.Parent.Accession), Is.EquivalentTo(new[] { "A", "B" }));
        Assert.That(index["EEEEER"].Select(p => p.Parent.Accession), Is.EquivalentTo(new[] { "A" }));
        Assert.That(index.Values.SelectMany(v => v).Count(p => p.Parent.IsDecoy), Is.GreaterThan(0), "decoy proteins are digested too");
    }

    /// <summary>
    /// Parsimony as in DDA: the shared peptide is explained by A, which has its own unique peptide, so B (seen only
    /// through the shared peptide) is not reported.
    /// </summary>
    [Test]
    public void ParsimonyExplainsASharedPeptideByTheProteinWithEvidenceOfItsOwn()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("LLLLLK", 0.001, index: 0),
            Peptide("EEEEER", 0.001, index: 1),
            Peptide("VVVVVK", 0.001, index: 2),
            Peptide("TTTTTR", 0.001, index: 3),
        };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        var names = result.ProteinGroups.Where(g => !g.IsDecoy).Select(g => g.ProteinGroupName).ToList();
        Assert.That(names, Is.EquivalentTo(new[] { "A", "C" }));
        var a = result.ProteinGroups.Single(g => g.ProteinGroupName == "A");
        Assert.That(a.AllPeptides.Select(p => p.FullSequence), Is.EquivalentTo(new[] { "LLLLLK", "EEEEER" }));
        // Unique among the groups reported: B is dropped, so the shared LLLLLK is A's alone (mzLib #1463, MetaMorpheus 1.1.13)
        Assert.That(a.UniquePeptides.Select(p => p.FullSequence), Is.EquivalentTo(new[] { "LLLLLK", "EEEEER" }));
        Assert.That(result.UnmappedPeptides, Is.EqualTo(0));
    }

    [Test]
    public void ADecoyPeptideMakesADecoyProteinGroup()
    {
        string decoy = DecoySequenceOf("C");
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("EEEEER", 0.001, score: 9, index: 0),
            Peptide(decoy, 0.002, score: 8, isDecoy: true, index: 1),
        };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        Assert.That(result.ProteinGroups.Single(g => g.IsDecoy).Proteins.Single().Accession, Is.EqualTo("DECOY_C"));
        Assert.That(result.ProteinGroups.Single(g => !g.IsDecoy).ProteinGroupName, Is.EqualTo("A"));
    }

    /// <summary>Only peptides that pass the q-value threshold count as evidence, exactly as DDA PSMs do.</summary>
    [Test]
    public void PeptidesAboveTheQValueThresholdAreNotEvidence()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("EEEEER", 0.001, index: 0),
            Peptide("VVVVVK", 0.5, index: 1),
        };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        Assert.That(result.ProteinGroups.Select(g => g.ProteinGroupName), Is.EquivalentTo(new[] { "A" }));
    }

    /// <summary>
    /// The DIA classifier's scores can be negative, and MetaMorpheus drops a spectral match scored at or below zero.
    /// A single order-preserving shift keeps every peptide, and the protein's best peptide is still the best scored.
    /// </summary>
    [Test]
    public void NegativeScoresAreShiftedNotDropped()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("EEEEER", 0.001, score: -1.5, index: 0),
            Peptide("LLLLLK", 0.001, score: -3.0, index: 1),
            Peptide("VVVVVK", 0.001, score: -7.25, index: 2),
        };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        Assert.That(result.SpectralMatches.Select(m => m.Score), Is.All.GreaterThan(0));
        Assert.That(result.SpectralMatches.OrderByDescending(m => m.Score).Select(m => m.FullSequence),
            Is.EqualTo(new[] { "EEEEER", "LLLLLK", "VVVVVK" }));
        Assert.That(result.ProteinGroups.Select(g => g.ProteinGroupName), Is.EquivalentTo(new[] { "A", "C" }));
    }

    [Test]
    public void EachSpectralMatchCarriesItsPeptideQValueAndPointsBackToItsPrecursor()
    {
        var peptides = new List<DiaPeptideMatch> { Peptide("EEEEER", 0.004, index: 41) };

        var match = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []).SpectralMatches.Single();

        Assert.That(match.FdrInfo.QValue, Is.EqualTo(0.004));
        Assert.That(match.FdrInfo.QValueNotch, Is.EqualTo(0.004));
        Assert.That(match.PeptideFdrInfo.QValue, Is.EqualTo(0.004));
        Assert.That(match.ScanNumber, Is.EqualTo(42), "one-based precursor index");
        Assert.That(match.FullFilePath, Is.EqualTo("run.raw"));
    }

    /// <summary>A peptide the database cannot yield is reported in a count, not silently lost.</summary>
    [Test]
    public void APeptideNoProteinYieldsIsCounted()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("EEEEER", 0.001, index: 0),
            Peptide("WWWWWK", 0.001, index: 1),
        };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        Assert.That(result.UnmappedPeptides, Is.EqualTo(1));
        Assert.That(result.SpectralMatches.Select(m => m.FullSequence), Is.EqualTo(new[] { "EEEEER" }));
    }

    /// <summary>
    /// A decoy peptide must land on a decoy protein, a target on a target. A sequence found only on the other side is
    /// unmapped rather than given a label it does not have.
    /// </summary>
    [Test]
    public void APeptideNeverLandsOnAProteinOfTheOtherKind()
    {
        var peptides = new List<DiaPeptideMatch> { Peptide("EEEEER", 0.001, isDecoy: true) };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        Assert.That(result.UnmappedPeptides, Is.EqualTo(1));
        Assert.That(result.ProteinGroups, Is.Empty);
    }

    // Peptide-level decoys, as a library made by reversing peptides has them: each pairs with one target sequence.
    private static readonly Dictionary<string, string> TargetOfDecoy = new()
    {
        ["KLLLLL"] = "LLLLLK",
        ["REEEEE"] = "EEEEER",
        ["KVVVVV"] = "VVVVVK",
        ["KWWWWW"] = "WWWWWK",
        ["RGGGGG"] = "GGGGGR",
    };

    private static DiaProteinInferenceResult InferWithInheritedDecoys(List<DiaPeptideMatch> peptides) =>
        DiaProteinInference.Infer(peptides, Targets, Parameters, "run.raw", [], [],
            decoy => TargetOfDecoy.GetValueOrDefault(decoy));

    /// <summary>
    /// A library whose decoys are reversed peptides has no decoy proteins to digest. Each decoy then takes its target's
    /// proteins, as decoys: the decoy of a shared peptide is shared by the same decoy proteins, so the decoy side has
    /// the target side's protein structure and protein FDR stays fair.
    /// </summary>
    [Test]
    public void ADecoyPeptideInheritsItsTargetsProteinsAsDecoys()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("KLLLLL", 0.001, isDecoy: true, index: 0),
            Peptide("REEEEE", 0.001, isDecoy: true, index: 1),
            Peptide("RGGGGG", 0.001, isDecoy: true, index: 2),
            Peptide("KVVVVV", 0.001, isDecoy: true, index: 3),
        };

        var result = InferWithInheritedDecoys(peptides);

        Assert.That(result.ProteinGroups.Select(g => g.ProteinGroupName), Is.EquivalentTo(new[] { "DECOY_A", "DECOY_B", "DECOY_C" }));
        foreach (string name in new[] { "DECOY_A", "DECOY_B" })
        {
            var group = result.ProteinGroups.Single(g => g.ProteinGroupName == name);
            Assert.That(group.AllPeptides.Select(p => p.FullSequence), Does.Contain("KLLLLL"), $"the shared decoy is in {name}");
            Assert.That(group.UniquePeptides.Select(p => p.FullSequence), Does.Not.Contain("KLLLLL"));
        }
        Assert.That(result.ProteinGroups, Is.All.Matches<ProteinGroup>(g => g.IsDecoy));
    }

    /// <summary>
    /// For decoys made from decoy proteins the database does not hold (permuted ones, say), the library's recorded
    /// accessions place them: each decoy goes to every decoy protein its entry names, and only there. Here the database
    /// has targets only, so digestion could place none of them.
    /// </summary>
    [Test]
    public void ADecoyPeptideMapsToTheDecoyProteinsItsLibraryEntryNames()
    {
        // mzLib's permuted decoys are named DECOY_<target>_fold<n> (EntrapmentAccession.Format)
        var accessions = new Dictionary<string, string[]>
        {
            ["ELEELK"] = ["DECOY_A_fold1", "DECOY_B_fold1"],
            ["QQQQQR"] = ["DECOY_A_fold1"],
            ["SSSSSK"] = ["DECOY_B_fold1"],
            ["NNNNNR"] = ["DECOY_C_fold1"],
        };
        var peptides = accessions.Keys.Select((s, i) => Peptide(s, 0.001, isDecoy: true, index: i)).ToList();

        var result = DiaProteinInference.Infer(peptides, Targets, Parameters, "run.raw", [], [],
            decoyAccessionsOf: s => accessions.GetValueOrDefault(s));

        Assert.That(result.UnmappedPeptides, Is.Zero);
        Assert.That(result.ProteinGroups, Is.All.Matches<ProteinGroup>(g => g.IsDecoy));
        Assert.That(result.ProteinGroups.Select(g => g.ProteinGroupName), Is.EquivalentTo(new[] { "DECOY_A_fold1", "DECOY_B_fold1", "DECOY_C_fold1" }));
        var a = result.ProteinGroups.Single(g => g.ProteinGroupName == "DECOY_A_fold1");
        Assert.That(a.AllPeptides.Select(p => p.FullSequence), Is.EquivalentTo(new[] { "ELEELK", "QQQQQR" }));
        Assert.That(a.AllPeptides.Select(p => p.Parent.BaseSequence.Substring(p.OneBasedStartResidue - 1, p.Length)), Is.EquivalentTo(new[] { "ELEELK", "QQQQQR" }),
            "each decoy sits in its decoy protein's sequence where its positions say");
        var b = result.ProteinGroups.Single(g => g.ProteinGroupName == "DECOY_B_fold1");
        Assert.That(b.AllPeptides.Select(p => p.FullSequence), Is.EquivalentTo(new[] { "ELEELK", "SSSSSK" }));
    }

    /// <summary>
    /// Decoy accessions not prefixed DECOY_ (A_DECOY, B_DECOY) are prefixed, because the picked protein FDR reads the text
    /// before the first '_' as the decoy identifier: unprefixed, every group would strip to "DECOY" and all but one be dropped.
    /// </summary>
    [Test]
    public void UnprefixedDecoyAccessionsArePrefixedSoPickedFdrKeepsThemApart()
    {
        var accessions = new Dictionary<string, string[]> { ["QQQQQR"] = ["A_DECOY"], ["SSSSSK"] = ["B_DECOY"], ["NNNNNR"] = ["C_DECOY"] };
        var peptides = accessions.Keys.Select((s, i) => Peptide(s, 0.001, isDecoy: true, index: i)).ToList();

        var result = DiaProteinInference.Infer(peptides, Targets, Parameters, "run.raw", [], [],
            decoyAccessionsOf: s => accessions.GetValueOrDefault(s));

        Assert.That(result.ProteinGroups.Select(g => g.ProteinGroupName), Is.EquivalentTo(new[] { "DECOY_A_DECOY", "DECOY_B_DECOY", "DECOY_C_DECOY" }));
    }

    /// <summary>A decoy whose library entry names no protein is not evidence, and is counted, as an unmapped target is.</summary>
    [Test]
    public void ADecoyItsLibraryNamesNoProteinForIsCounted()
    {
        var peptides = new List<DiaPeptideMatch> { Peptide("ELEELK", 0.001, isDecoy: true, index: 0), Peptide("QQQQQR", 0.001, isDecoy: true, index: 1) };

        var result = DiaProteinInference.Infer(peptides, Targets, Parameters, "run.raw", [], [],
            decoyAccessionsOf: s => s == "ELEELK" ? ["A_DECOY"] : null);

        Assert.That(result.UnmappedPeptides, Is.EqualTo(1));
        Assert.That(result.SpectralMatches.Select(m => m.FullSequence), Is.EquivalentTo(new[] { "ELEELK" }));
    }

    /// <summary>
    /// Picked protein FDR (PLAN D4) on a hand-built pair: a protein and its decoy compete, and only the better of the two
    /// counts. With the target better, its decoy is picked away and the target passes; with the decoy better, the decoy is
    /// kept. MetaMorpheus then adds the losing target back with its classic q-value; that value is not asserted, because
    /// the classic pass sorts on BestPeptideQValue and BestPeptideScore before they are set (ProteinScoringAndFdrEngine
    /// assigns them only in the picked loop), so it ranks groups in input order and here gives the losing target 0.
    /// </summary>
    [TestCase(false)]
    [TestCase(true)]
    public void PickedProteinFdrKeepsTheBetterOfATargetAndItsDecoy(bool decoyBetter)
    {
        string decoy = DecoySequenceOf("A");
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("EEEEER", decoyBetter ? 0.002 : 0.001, score: decoyBetter ? 8 : 10, index: 0),
            Peptide(decoy, decoyBetter ? 0.001 : 0.002, score: decoyBetter ? 10 : 8, isDecoy: true, index: 1),
        };

        var result = DiaProteinInference.Infer(peptides, TargetsAndDecoys(), Parameters, "run.raw", [], []);

        var names = result.ProteinGroups.Select(g => g.ProteinGroupName).ToList();
        if (decoyBetter)
        {
            Assert.That(names, Does.Contain("DECOY_A"), "the better decoy is kept");
        }
        else
        {
            Assert.That(names, Does.Not.Contain("DECOY_A"), "the worse decoy is picked away");
            Assert.That(result.ProteinGroups.Single(g => g.ProteinGroupName == "A").QValue, Is.Zero);
        }
    }

    [Test]
    public void DecoysAreMappedOneWayOnly()
    {
        var e = Assert.Throws<ArgumentException>(() => DiaProteinInference.Infer([], Targets, Parameters, "run.raw", [], [],
            targetSequenceOfDecoy: _ => null, decoyAccessionsOf: _ => null));
        Assert.That(e!.Message, Does.Contain("decoy"));
    }

    /// <summary>One decoy protein per target protein, whichever of its decoy peptides is seen; and with no other evidence, the shared decoy is explained by DECOY_A, as its target is by A.</summary>
    [Test]
    public void DecoysOfOneProteinShareOneDecoyProtein()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("KLLLLL", 0.001, isDecoy: true, index: 0),
            Peptide("REEEEE", 0.001, isDecoy: true, index: 1),
        };

        var result = InferWithInheritedDecoys(peptides);

        var a = result.ProteinGroups.Single(g => g.ProteinGroupName == "DECOY_A");
        Assert.That(a.Proteins, Has.Count.EqualTo(1));
        Assert.That(a.AllPeptides.Select(p => p.FullSequence), Is.EquivalentTo(new[] { "KLLLLL", "REEEEE" }));
    }

    [Test]
    public void ADecoyWhoseTargetNoProteinYieldsIsCounted()
    {
        var peptides = new List<DiaPeptideMatch>
        {
            Peptide("KWWWWW", 0.001, isDecoy: true, index: 0),
            Peptide("KNOPAIR", 0.001, isDecoy: true, index: 1),
            Peptide("EEEEER", 0.001, index: 2),
        };

        var result = InferWithInheritedDecoys(peptides);

        Assert.That(result.UnmappedPeptides, Is.EqualTo(2));
        Assert.That(result.ProteinGroups.Select(g => g.ProteinGroupName), Is.EquivalentTo(new[] { "A" }));
    }

    /// <summary>
    /// The DIA search has no PEP. A PEP q-value threshold would switch the protein engine to filtering on PEP, which
    /// every DIA peptide fails, and report no proteins at all. Refuse it instead.
    /// </summary>
    [Test]
    public void APepQValueThresholdIsRefused()
    {
        var pepFiltered = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5), pepQValueThreshold: 0.005);

        Assert.Throws<ArgumentException>(() =>
            DiaProteinInference.Infer([Peptide("EEEEER", 0.001)], TargetsAndDecoys(), pepFiltered, "run.raw", [], []));
    }
}
