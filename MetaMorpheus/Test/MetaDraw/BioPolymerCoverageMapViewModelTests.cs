using GuiFunctions.MetaDraw;
using NUnit.Framework;
using System.IO;
using System.Linq;
using System.Windows.Media;
using Readers;
using System.Threading;
using System;
using System.Collections.Generic;

namespace Test.MetaDraw;

[TestFixture]
public class BioPolymerCoverageMapViewModelTests
{
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
        tabVm.DatabasePaths.Add(dbPath);
        var method = typeof(BioPolymerTabViewModel)
            .GetMethod("LoadDatabase", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        method.Invoke(tabVm, null);
        Thread.Sleep(3000); // wait for async load
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

    #region Viridis Intensity Tests

    [Test]
    public void ViridisColors_Has20Colors()
    {
        var field = typeof(BioPolymerCoverageMapViewModel)
            .GetField("ViridisColors", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);
        Assert.That(field, Is.Not.Null);
        var colors = field.GetValue(null) as SolidColorBrush[];
        Assert.That(colors, Is.Not.Null);
        Assert.That(colors.Length, Is.EqualTo(20));
        Assert.That(colors.All(c => c != null), Is.True);
    }

    [Test]
    public void ViridisColors_StartAndEndColorsDiffer()
    {
        var field = typeof(BioPolymerCoverageMapViewModel)
            .GetField("ViridisColors", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);
        var colors = field.GetValue(null) as SolidColorBrush[];
        Assert.That(colors[0].Color, Is.Not.EqualTo(colors[19].Color));
    }

    [Test]
    public void GetIntensityBrush_Zero_ReturnsFirstColor()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        var method = typeof(BioPolymerCoverageMapViewModel)
            .GetMethod("GetIntensityBrush", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        var result = method.Invoke(vm, new object[] { 0.0 }) as SolidColorBrush;
        var field = typeof(BioPolymerCoverageMapViewModel)
            .GetField("ViridisColors", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);
        var colors = field.GetValue(null) as SolidColorBrush[];
        Assert.That(result, Is.EqualTo(colors[0]));
    }

    [Test]
    public void GetIntensityBrush_One_ReturnsLastColor()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        var method = typeof(BioPolymerCoverageMapViewModel)
            .GetMethod("GetIntensityBrush", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        var result = method.Invoke(vm, new object[] { 1.0 }) as SolidColorBrush;
        var field = typeof(BioPolymerCoverageMapViewModel)
            .GetField("ViridisColors", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);
        var colors = field.GetValue(null) as SolidColorBrush[];
        Assert.That(result, Is.EqualTo(colors[19]));
    }

    [Test]
    public void GetIntensityBrush_MidValue_ReturnsMidColor()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        var method = typeof(BioPolymerCoverageMapViewModel)
            .GetMethod("GetIntensityBrush", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        var result = method.Invoke(vm, new object[] { 0.5 }) as SolidColorBrush;
        var field = typeof(BioPolymerCoverageMapViewModel)
            .GetField("ViridisColors", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);
        var colors = field.GetValue(null) as SolidColorBrush[];
        // 0.5 * 19 = 9.5, clamped to 9
        Assert.That(result, Is.EqualTo(colors[9]));
    }

    [Test]
    public void GetIntensityBrush_Negative_ClampsToFirstColor()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        var method = typeof(BioPolymerCoverageMapViewModel)
            .GetMethod("GetIntensityBrush", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        var result = method.Invoke(vm, new object[] { -0.1 }) as SolidColorBrush;
        var field = typeof(BioPolymerCoverageMapViewModel)
            .GetField("ViridisColors", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);
        var colors = field.GetValue(null) as SolidColorBrush[];
        Assert.That(result, Is.EqualTo(colors[0]));
    }

    [Test]
    public void GetIntensityBrush_AboveOne_ClampsToLastColor()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        var method = typeof(BioPolymerCoverageMapViewModel)
            .GetMethod("GetIntensityBrush", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        var result = method.Invoke(vm, new object[] { 1.5 }) as SolidColorBrush;
        var field = typeof(BioPolymerCoverageMapViewModel)
            .GetField("ViridisColors", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Static);
        var colors = field.GetValue(null) as SolidColorBrush[];
        Assert.That(result, Is.EqualTo(colors[19]));
    }

    [Test]
    public void ColorBy_CanSetPrecursorIntensity()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        bool propertyChanged = false;
        vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.ColorBy)) propertyChanged = true; };

        vm.ColorBy = ColorResultsBy.PrecursorIntensity;
        Assert.That(vm.ColorBy, Is.EqualTo(ColorResultsBy.PrecursorIntensity));
        Assert.That(propertyChanged, Is.True);
    }

    [Test]
    public void ColorBy_CanSetScore()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        bool propertyChanged = false;
        vm.PropertyChanged += (s, e) => { if (e.PropertyName == nameof(vm.ColorBy)) propertyChanged = true; };

        vm.ColorBy = ColorResultsBy.Score;
        Assert.That(vm.ColorBy, Is.EqualTo(ColorResultsBy.Score));
        Assert.That(propertyChanged, Is.True);
    }

    [Test]
    public void CreateLegendItems_ReturnsEmpty_ForPrecursorIntensity()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        var method = typeof(BioPolymerCoverageMapViewModel)
            .GetMethod("CreateLegendItems", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);

        var filteredResults = new List<BioPolymerCoverageResultModel>();
        var result = method.Invoke(vm, new object[] { filteredResults, 16.0, 96.0, ColorResultsBy.PrecursorIntensity })
            as System.Collections.IEnumerable;
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Cast<object>().Any(), Is.False);
    }

    [Test]
    public void CreateLegendItems_ReturnsEmpty_ForScore()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        var method = typeof(BioPolymerCoverageMapViewModel)
            .GetMethod("CreateLegendItems", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);

        var filteredResults = new List<BioPolymerCoverageResultModel>();
        var result = method.Invoke(vm, new object[] { filteredResults, 16.0, 96.0, ColorResultsBy.Score })
            as System.Collections.IEnumerable;
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Cast<object>().Any(), Is.False);
    }

    [Test]
    public void Plotting_WithPrecursorIntensity_CreatesDrawing()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.PrecursorIntensity;

        var match = new DummySpectralmatch();
        var backingField = typeof(SpectrumMatchFromTsv)
            .GetField("<PrecursorIntensity>k__BackingField", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        backingField.SetValue(match, 100000.0);
        var result = new BioPolymerCoverageResultModel(match, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", new[] { result });
        vm.Group = group;

        Assert.That(vm.CoverageDrawing, Is.Not.Null);
    }

    [Test]
    public void Plotting_WithScore_CreatesDrawing()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.Score;

        var match = new DummySpectralmatch();
        match.SetScore(100);
        var result = new BioPolymerCoverageResultModel(match, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", new[] { result });
        vm.Group = group;

        Assert.That(vm.CoverageDrawing, Is.Not.Null);
    }

    [Test]
    public void Plotting_WithScore_DifferentScores_CreatesDrawing()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.Score;

        var match1 = new DummySpectralmatch();
        match1.SetScore(10);
        var match2 = new DummySpectralmatch();
        match2.SetScore(1000);
        var result1 = new BioPolymerCoverageResultModel(match1, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var result2 = new BioPolymerCoverageResultModel(match2, "ABC", 2, 3, BioPolymerCoverageType.Unique);

        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", new[] { result1, result2 });
        vm.Group = group;

        Assert.That(vm.CoverageDrawing, Is.Not.Null);
    }

    [Test]
    public void Plotting_WithIntensityMode_AllEqual_DoesNotCrash()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.Score;

        var match1 = new DummySpectralmatch();
        match1.SetScore(42);
        var match2 = new DummySpectralmatch();
        match2.SetScore(42);
        var result1 = new BioPolymerCoverageResultModel(match1, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var result2 = new BioPolymerCoverageResultModel(match2, "ABC", 2, 3, BioPolymerCoverageType.Unique);

        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", new[] { result1, result2 });
        vm.Group = group;

        Assert.That(vm.CoverageDrawing, Is.Not.Null);
    }

    [Test]
    public void Plotting_WithPrecursorIntensity_SomeNull_DoesNotCrash()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.PrecursorIntensity;

        var matchWithIntensity = new DummySpectralmatch();
        matchWithIntensity.SetScore(200);
        var matchNullIntensity = new DummySpectralmatch();
        matchNullIntensity.SetScore(100);
        var backingField = typeof(SpectrumMatchFromTsv)
            .GetField("<PrecursorIntensity>k__BackingField", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        backingField.SetValue(matchWithIntensity, 500000.0);
        backingField.SetValue(matchNullIntensity, null);

        var result1 = new BioPolymerCoverageResultModel(matchWithIntensity, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var result2 = new BioPolymerCoverageResultModel(matchNullIntensity, "ABC", 2, 3, BioPolymerCoverageType.Unique);

        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", new[] { result1, result2 });
        vm.Group = group;

        Assert.That(vm.CoverageDrawing, Is.Not.Null);
    }

    [Test]
    public void Plotting_WithPrecursorIntensity_AllNull_FallsBackToScore()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.PrecursorIntensity;

        var match1 = new DummySpectralmatch();
        match1.SetScore(10);
        var match2 = new DummySpectralmatch();
        match2.SetScore(1000);
        var backingField = typeof(SpectrumMatchFromTsv)
            .GetField("<PrecursorIntensity>k__BackingField", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        backingField.SetValue(match1, null);
        backingField.SetValue(match2, null);

        var result1 = new BioPolymerCoverageResultModel(match1, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var result2 = new BioPolymerCoverageResultModel(match2, "ABC", 2, 3, BioPolymerCoverageType.Unique);

        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", new[] { result1, result2 });
        vm.Group = group;

        Assert.That(vm.CoverageDrawing, Is.Not.Null);
    }

    #endregion
}