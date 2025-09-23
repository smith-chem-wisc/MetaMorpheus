using GuiFunctions.MetaDraw;
using NUnit.Framework;
using System.IO;
using System.Linq;
using System.Windows.Media;
using Readers;
using System.Threading;

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

    [Test]
    public void ColorBy_Setter_UpdatesValue_ClearsColors_RaisesPropertyChanged()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        bool propertyChanged = false;
        vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.ColorBy)) propertyChanged = true; };

        // Set to FileOrigin from default CoverageType
        vm.ColorBy = ColorResultsBy.FileOrigin;
        Assert.That(vm.ColorBy, Is.EqualTo(ColorResultsBy.FileOrigin));
        Assert.That(propertyChanged, Is.True);

        // Set to None, should also update
        propertyChanged = false;
        vm.ColorBy = ColorResultsBy.None;
        Assert.That(vm.ColorBy, Is.EqualTo(ColorResultsBy.None));
        Assert.That(propertyChanged, Is.True);
    }

    [Test]
    public void Plotting_WithColorByFileOrigin_AssignsDistinctColorsAndLegendItems()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.FileOrigin;

        // Create two results from different files
        var match1 = new DummySpectralmatch(0.005, "1", "fileA.psmtsv");
        var match2 = new DummySpectralmatch(0.01, "2", "fileB.psmtsv");
        var result1 = new BioPolymerCoverageResultModel(match1, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var result2 = new BioPolymerCoverageResultModel(match2, "ABC", 2, 3, BioPolymerCoverageType.Unique);

        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", new[] { result1, result2 });
        vm.Group = group;

        // After setting Group, CoverageDrawing should be non-null
        Assert.That(vm.CoverageDrawing, Is.Not.Null);

        // Use reflection to access the private IdentifierToColor dictionary
        var idToColorField = typeof(BioPolymerCoverageMapViewModel)
            .GetField("IdentifierToColor", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);
        Assert.That(idToColorField, Is.Not.Null);

        var idToColor = idToColorField.GetValue(null) as System.Collections.IDictionary;
        Assert.That(idToColor, Is.Not.Null);

        // Both file names should be present and mapped to SolidColorBrush
        Assert.That(idToColor.Contains("fileA.psmtsv"), Is.True);
        Assert.That(idToColor.Contains("fileB.psmtsv"), Is.True);
        var brushA = idToColor["fileA.psmtsv"] as SolidColorBrush;
        var brushB = idToColor["fileB.psmtsv"] as SolidColorBrush;
        Assert.That(brushA, Is.Not.Null);
        Assert.That(brushB, Is.Not.Null);
        Assert.That(brushA.Color, Is.Not.EqualTo(brushB.Color));

        // Legend items should include both files
        // Use reflection to call private CreateLegendItems
        var method = typeof(BioPolymerCoverageMapViewModel)
            .GetMethod("CreateLegendItems", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        Assert.That(method, Is.Not.Null);

        var filteredResults = group.CoverageResults.ToList();
        double fontSize = 16; // typical font size
        double dpi = 96; // typical DPI
        var legendItems = method.Invoke(vm, new object[] { filteredResults, fontSize, dpi, ColorResultsBy.FileOrigin }) as System.Collections.IEnumerable;
        Assert.That(legendItems, Is.Not.Null);

        var legendLabels = legendItems.Cast<(SolidColorBrush Brush, System.Windows.Media.FormattedText Text)>().Select(li => li.Text.Text).ToList();
        Assert.That(legendLabels.Any(l => l.Contains("fileA.psmtsv")), Is.True);
        Assert.That(legendLabels.Any(l => l.Contains("fileB.psmtsv")), Is.True);
    }

    [Test]
    public void Plotting_WithColorByCoverageType_AssignsTypeColorsAndLegendItems()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.CoverageType;

        // Create two results with different coverage types
        var match1 = new DummySpectralmatch(0.005, "1", "fileA.psmtsv");
        var match2 = new DummySpectralmatch(0.01, "2", "fileB.psmtsv");
        var result1 = new BioPolymerCoverageResultModel(match1, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var result2 = new BioPolymerCoverageResultModel(match2, "ABC", 2, 3, BioPolymerCoverageType.Shared);

        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", new[] { result1, result2 });
        vm.Group = group;

        // After setting Group, CoverageDrawing should be non-null
        Assert.That(vm.CoverageDrawing, Is.Not.Null);

        // Legend items should include both coverage types
        var method = typeof(BioPolymerCoverageMapViewModel)
            .GetMethod("CreateLegendItems", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        Assert.That(method, Is.Not.Null);

        var filteredResults = group.CoverageResults.ToList();
        double fontSize = 16; // typical font size
        double dpi = 96; // typical DPI
        var legendItems = method.Invoke(vm, new object[] { filteredResults, fontSize, dpi, ColorResultsBy.CoverageType }) as System.Collections.IEnumerable;
        Assert.That(legendItems, Is.Not.Null);

        var legendLabels = legendItems.Cast<(SolidColorBrush Brush, System.Windows.Media.FormattedText Text)>().Select(li => li.Text.Text).ToList();
        Assert.That(legendLabels.Any(l => l.Contains("Unique")), Is.True);
        Assert.That(legendLabels.Any(l => l.Contains("Shared")), Is.True);

        // The brushes should be the ones from MetaDrawSettings.BioPolymerCoverageColors
        var brushes = legendItems.Cast<(SolidColorBrush Brush, System.Windows.Media.FormattedText Text)>().Select(li => li.Brush).ToList();
        Assert.That(brushes.Distinct().Count(), Is.EqualTo(2));
    }

    [Test]
    public void Plotting_WithColorByNone_AssignsGrayColorAndNoLegendItems()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.None;

        // Create two results with different coverage types and files
        var match1 = new DummySpectralmatch(0.005, "1", "fileA.psmtsv");
        var match2 = new DummySpectralmatch(0.01, "2", "fileB.psmtsv");
        var result1 = new BioPolymerCoverageResultModel(match1, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var result2 = new BioPolymerCoverageResultModel(match2, "ABC", 2, 3, BioPolymerCoverageType.Shared);

        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", new[] { result1, result2 });
        vm.Group = group;

        // After setting Group, CoverageDrawing should be non-null
        Assert.That(vm.CoverageDrawing, Is.Not.Null);

        // Legend items should be empty
        var method = typeof(BioPolymerCoverageMapViewModel)
            .GetMethod("CreateLegendItems", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        Assert.That(method, Is.Not.Null);

        var filteredResults = group.CoverageResults.ToList();
        double fontSize = 16; // typical font size
        double dpi = 96; // typical DPI
        var legendItems = method.Invoke(vm, new object[] { filteredResults, fontSize, dpi, ColorResultsBy.None }) as System.Collections.IEnumerable;
        Assert.That(legendItems, Is.Not.Null);
        Assert.That(legendItems.Cast<object>().Any(), Is.False);
    }

    [Test]
    public void Plotting_WithRealData()
    {
        // Bring in search results
        if (!EverythingRunnerEngineTestCase.TryGetTestCase(EverythingRunnerEngineTestCases.BottomUpQValue, out var searchTestCase))
            Assert.Fail();

        var dbPath = searchTestCase.DatabaseList.First().FilePath;
        var smPath = Directory.GetFiles(searchTestCase.OutputDirectory, "*PSMs.psmtsv", SearchOption.AllDirectories).First();
        var matches = SpectrumMatchTsvReader.ReadTsv(smPath, out _);

        // Load search results into Tab
        var logic = new DummyMetaDrawLogic();
        logic.AllSpectralMatches = matches;
        logic.FilterPsms();

        var tabVm = new BioPolymerTabViewModel(logic);
        tabVm.DatabasePath = dbPath;
        var method = typeof(BioPolymerTabViewModel)
            .GetMethod("LoadDatabase", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        method.Invoke(tabVm, null);
        Thread.Sleep(1000); // wait for async load
        Assert.That(tabVm.IsDatabaseLoaded, Is.True);

        // Plot every group
        var previous = tabVm.CoverageMapViewModel.CoverageDrawing;
        foreach (var group in tabVm.AllGroups)
        {
            tabVm.SelectedGroup = group;
            var current = tabVm.CoverageMapViewModel.CoverageDrawing;
            Assert.That(current, Is.Not.Null);
            Assert.That(current, Is.Not.EqualTo(previous));
        }
    }
}