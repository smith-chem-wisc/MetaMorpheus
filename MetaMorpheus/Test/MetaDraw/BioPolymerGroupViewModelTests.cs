using GuiFunctions;
using GuiFunctions.MetaDraw;
using NUnit.Framework;

namespace Test.MetaDraw;

[TestFixture]
public class BioPolymerGroupViewModelTests
{
    [Test]
    public void Constructor_InitializesProperties()
    {
        var match = new DummySpectralmatch(0.01, "1", "f.psmtsv");
        var results = new[]
        {
            new BioPolymerCoverageResultModel(match, "ABC", 1, 2, BioPolymerCoverageType.Unique),
            new BioPolymerCoverageResultModel(match, "ABC", 2, 3, BioPolymerCoverageType.Shared)
        };
        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", results);

        Assert.That(group.Accession, Is.EqualTo("ACC"));
        Assert.That(group.ProteinName, Is.EqualTo("Prot"));
        Assert.That(group.Sequence, Is.EqualTo("ABC"));
        Assert.That(group.Length, Is.EqualTo(3));
        Assert.That(group.DatabaseName, Is.EqualTo(string.Empty));
        Assert.That(group.CoverageResults.Count, Is.EqualTo(2));
    }

    [Test]
    public void UpdatePropertiesAfterFilter_ComputesCoverage()
    {
        var match = new DummySpectralmatch(0.01, "1", "f.psmtsv");
        var results = new[]
        {
            new BioPolymerCoverageResultModel(match, "ABC", 1, 2, BioPolymerCoverageType.Unique),
            new BioPolymerCoverageResultModel(match, "ABC", 2, 3, BioPolymerCoverageType.UniqueMissedCleavage),
            new BioPolymerCoverageResultModel(match, "ABC", 1, 3, BioPolymerCoverageType.Shared)
        };
        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", results, "test.fasta");
        group.UpdatePropertiesAfterFilter();

        Assert.That(group.GroupCount, Is.EqualTo(3));
        Assert.That(group.DatabaseName, Is.EqualTo("test.fasta"));
        Assert.That(group.UniqueSequenceCoverage, Is.GreaterThan(0));
        Assert.That(group.MaximumSequenceCoverage, Is.EqualTo(1.0));
    }

    [Test]
    public void UpdatePropertiesAfterFilter_EmptyResults()
    {
        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", new BioPolymerCoverageResultModel[0]);
        group.UpdatePropertiesAfterFilter();
        Assert.That(group.GroupCount, Is.EqualTo(0));
        Assert.That(group.UniqueSequenceCoverage, Is.EqualTo(0));
        Assert.That(group.MaximumSequenceCoverage, Is.EqualTo(0));
    }
}