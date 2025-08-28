using GuiFunctions;
using NUnit.Framework;

namespace Test.MetaDraw;

[TestFixture]
public class BioPolymerCoverageMapViewModelTests
{
    [Test]
    public void CoverageDrawing_Null_WhenGroupIsNull()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.Group = null;
        Assert.That(vm.CoverageDrawing, Is.Null);
    }

    [Test]
    public void Group_Set_TriggersRedraw()
    {
        var match = new DummySpectralmatch(0.01, "1", "f.psmtsv");
        var results = new[]
        {
            new BioPolymerCoverageResultModel(match, "ABC", 1, 2, BioPolymerCoverageType.Unique)
        };
        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", results);
        var vm = new BioPolymerCoverageMapViewModel();
        vm.Group = group;
        Assert.That(vm.Group, Is.EqualTo(group));
        Assert.That(vm.CoverageDrawing, Is.Not.Null);
    }

    [Test]
    public void UpdateLettersPerRow_ChangesLettersPerRow()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        int before = vm.LettersPerRow;
        vm.UpdateLettersPerRow(2000);
        Assert.That(vm.LettersPerRow, Is.Not.EqualTo(before));
    }

    [Test]
    public void Redraw_WithNoResults_DoesNotThrow()
    {
        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", new BioPolymerCoverageResultModel[0]);
        var vm = new BioPolymerCoverageMapViewModel();
        vm.Group = group;
        Assert.That(vm.CoverageDrawing, Is.Not.Null);
    }
}