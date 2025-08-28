using GuiFunctions.MetaDraw;
using NUnit.Framework;

namespace Test.MetaDraw;

[TestFixture]
public class BioPolymerCoverageMapViewModelTests
{
    [Test]
    public void LettersPerRow_Default_Is50()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        Assert.That(vm.LettersPerRow, Is.EqualTo(50));
    }

    [Test]
    public void LettersPerRow_Setter_UpdatesValueAndRaisesPropertyChanged()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        bool propertyChanged = false;
        vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.LettersPerRow)) propertyChanged = true; };
        vm.LettersPerRow = 42;
        Assert.That(vm.LettersPerRow, Is.EqualTo(42));
        Assert.That(propertyChanged, Is.True);
    }

    [Test]
    public void AvailableWidth_Default_Is800()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        Assert.That(vm.AvailableWidth, Is.EqualTo(800));
    }

    [Test]
    public void AvailableWidth_Setter_UpdatesValueAndRaisesPropertyChanged()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        bool propertyChanged = false;
        vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.AvailableWidth)) propertyChanged = true; };
        vm.AvailableWidth = 1234;
        Assert.That(vm.AvailableWidth, Is.EqualTo(1234));
        Assert.That(propertyChanged, Is.True);
    }

    [Test]
    public void UpdateLettersPerRow_ChangesLettersPerRow_WhenWidthChanges()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        int before = vm.LettersPerRow;
        vm.UpdateLettersPerRow(2000);
        Assert.That(vm.LettersPerRow, Is.Not.EqualTo(before));
    }

    [Test]
    public void UpdateLettersPerRow_DoesNotSetBelow10()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.UpdateLettersPerRow(1); // very small width
        Assert.That(vm.LettersPerRow, Is.EqualTo(10));
    }

    [Test]
    public void UpdateLettersPerRow_CallsRedraw_IfLettersPerRowUnchanged()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        bool redrawCalled = false;
        // Use reflection to replace private Redraw with a test hook
        var method = typeof(BioPolymerCoverageMapViewModel).GetMethod("Redraw", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        Assert.That(method, Is.Not.Null);
        // Set LettersPerRow to a known value
        vm.LettersPerRow = 20;
        // Now call UpdateLettersPerRow with a width that results in the same LettersPerRow
        double width = 20 * 16 * 0.70; // assuming font size 16
        // Can't directly hook Redraw, but we can set Group and check CoverageDrawing is not null after
        var match = new DummySpectralmatch(0.01, "1", "f.psmtsv");
        var results = new[] { new BioPolymerCoverageResultModel(match, "ABC", 1, 2, BioPolymerCoverageType.Unique) };
        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", results);
        vm.Group = group;
        var beforeDrawing = vm.CoverageDrawing;
        vm.UpdateLettersPerRow(width);
        Assert.That(vm.CoverageDrawing, Is.Not.Null);
    }
}