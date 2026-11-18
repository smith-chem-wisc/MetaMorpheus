#nullable enable
using System;
using System.Diagnostics.CodeAnalysis;
using EngineLayer.DiaLibrarySearch;
using NUnit.Framework;

namespace Test.DiaLibrarySearch;

/// <summary>
/// Protein-level support: how many other peptides of a precursor's protein are confident. DIA-NN 2.3.2 re-scores its
/// precursors on protein information after its network (+4.5% at 1% q on PXD071205 A1 with each entrapment in its owner's
/// protein). A decoy carries its target's protein, so it gets the support a false hit in that protein would.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaProteinSupportTests
{
    private static DiaProteinSupport.Row Target(string accession, string sequence, bool confident) => new(accession, sequence, IsDecoy: false, confident);
    private static DiaProteinSupport.Row Decoy(string accession, string sequence) => new(accession, sequence, IsDecoy: true, Confident: false);

    [Test]
    public void SupportCountsTheOtherConfidentPeptidesOfAPrecursorsProteins()
    {
        DiaProteinSupport.Row[] rows =
        [
            Target("P1", "AAAK", confident: true),     // 0
            Target("P1", "AAAK", confident: true),     // 1: the same peptide at another charge
            Target("P1", "CCCK", confident: true),     // 2
            Target("P1", "DDDK", confident: false),    // 3
            Target("P2", "EEEK", confident: true),     // 4
            Target("P2|P1", "FFFK", confident: false), // 5: shared, read on its best-supported protein
            Decoy("DECOY_P1", "KAAA"),                 // 6: a decoy counts toward its target's protein
            Decoy("DECOY_P2", "KEEE"),                 // 7
            Target("Random_P1_f0", "AKAA", confident: false), // 8: an entrapment in its own protein
            Target("P3", "GGGK", confident: false),    // 9: a protein with nothing confident
            Target("", "HHHK", confident: true),       // 10: no accession
        ];
        double[] support = DiaProteinSupport.Measure(rows);
        // One peptide, however many charges; never the row's own peptide
        Assert.That(support, Is.EqualTo(new double[] { 1, 1, 1, 2, 0, 2, 2, 1, 0, 0, 0 }));

        // Each entrapment counted toward its owner's protein: support as a false hit in a present protein would get
        static string Owner(string accession) => accession.StartsWith("Random_", StringComparison.Ordinal) ? accession[7..accession.LastIndexOf("_f", StringComparison.Ordinal)] : accession;
        Assert.That(DiaProteinSupport.Measure(rows, Owner)[8], Is.EqualTo(2));
        Assert.That(DiaProteinSupport.Measure(rows, Owner)[..8], Is.EqualTo(support[..8]));
    }

    [Test]
    public void ProteinSupportIsOffByDefaultAndCalibrationLeavesItOff()
    {
        Assert.That(new DiaLibrarySearchParameters().ProteinSupport, Is.False);
        Assert.That(new DiaLibrarySearchParameters().ProteinKey, Is.Null);
        Assert.That(DiaIrtSelfCalibration.CalibrationParameters(new DiaLibrarySearchParameters(ProteinSupport: true)).ProteinSupport, Is.False);
        Assert.That(DiaNeighbourEvidence.ProteinSupportColumnName, Is.EqualTo("ProteinSupport"));
    }
}
