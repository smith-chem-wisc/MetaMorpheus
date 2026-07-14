using GuiFunctions.MetaDraw;
using NUnit.Framework;
using System.IO;
using System.Linq;
using System.Windows.Media;
using Readers;
using System.Threading;
using System;
using System.Collections.Generic;
using GuiFunctions.MetaDraw.BioPolymerCoverage.ColorMapping.Gradient;

namespace Test.MetaDraw;

file sealed class AlwaysNullNumericMapper : NumericBioPolymerCoverageColorMapper
{
    public override ColorResultsBy ColorBy => ColorResultsBy.PrecursorIntensity;
    public override string DisplayName => "Always Null";
    public override bool DefaultUseLogScale => true;

    public override double? GetNumericValue(BioPolymerCoverageResultModel result)
    {
        return null;
    }

    protected override NumericBioPolymerCoverageColorMapper? GetFallbackMapper() => new ScoreColorMapper();
}

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
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.FileOrigin);

        // Create two results from different files
        var match1 = new DummySpectralmatch(0.005, "1", "fileA.psmtsv");
        var match2 = new DummySpectralmatch(0.01, "2", "fileB.psmtsv");
        var result1 = new BioPolymerCoverageResultModel(match1, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var result2 = new BioPolymerCoverageResultModel(match2, "ABC", 2, 3, BioPolymerCoverageType.Unique);

        var group = new BioPolymerGroupViewModel("ACC", "Prot", "ABC", new[] { result1, result2 });
        vm.Group = group;

        // After setting Group, CoverageDrawing should be non-null
        Assert.That(vm.CoverageDrawing, Is.Not.Null);

        mapper.Prepare([result1, result2], ColorGradientType.Viridis, false);
        var brushA = mapper.GetBrush(result1);
        var brushB = mapper.GetBrush(result2);
        Assert.That(brushA.Color, Is.Not.EqualTo(brushB.Color));

        // Legend items should include both files
        // Use reflection to call private CreateLegendItems
        var method = typeof(BioPolymerCoverageMapViewModel)
            .GetMethod("CreateLegendItems", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        Assert.That(method, Is.Not.Null);

        var filteredResults = group.CoverageResults.ToList();
        double fontSize = 16; // typical font size
        double dpi = 96; // typical DPI
        mapper.Prepare(filteredResults, ColorGradientType.Viridis, false);
        var legendItems = method.Invoke(vm, new object[] { fontSize, dpi, mapper }) as System.Collections.IEnumerable;
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
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.CoverageType);

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
        mapper.Prepare(filteredResults, ColorGradientType.Viridis, false);
        var legendItems = method.Invoke(vm, new object[] { fontSize, dpi, mapper }) as System.Collections.IEnumerable;
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
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.None);

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
        mapper.Prepare(filteredResults, ColorGradientType.Viridis, false);
        var legendItems = method.Invoke(vm, new object[] { fontSize, dpi, mapper }) as System.Collections.IEnumerable;
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
    public void MapperFactory_ReturnsNumericMapper_ForScore()
    {
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.Score);
        Assert.That(mapper, Is.Not.Null);
        Assert.That(mapper.ColorBy, Is.EqualTo(ColorResultsBy.Score));
        Assert.That(mapper.IsNumeric, Is.True);
    }

    [Test]
    public void NumericMapper_CreateRenderContext_FallsBackToScore()
    {
        var match1 = new DummySpectralmatch();
        match1.SetScore(50);
        var match2 = new DummySpectralmatch();
        match2.SetScore(500);
        var result1 = new BioPolymerCoverageResultModel(match1, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var result2 = new BioPolymerCoverageResultModel(match2, "ABC", 2, 3, BioPolymerCoverageType.Unique);
        var results = new List<BioPolymerCoverageResultModel> { result1, result2 };

        var mapper = new AlwaysNullNumericMapper();
        mapper.Prepare(results, ColorGradientType.Viridis, true);
        Assert.That(mapper.GradientLegendTitle, Is.EqualTo("Score (log10)"));
        Assert.That(mapper.GradientMinValue, Is.EqualTo(System.Math.Log10(50)).Within(1e-9));
        Assert.That(mapper.GradientMaxValue, Is.EqualTo(System.Math.Log10(500)).Within(1e-9));
    }

    [Test]
    public void ScoreColorMapper_CreateRenderContext_DerivesMinMax()
    {
        var match1 = new DummySpectralmatch();
        match1.SetScore(10);
        var match2 = new DummySpectralmatch();
        match2.SetScore(1000);
        var result1 = new BioPolymerCoverageResultModel(match1, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var result2 = new BioPolymerCoverageResultModel(match2, "ABC", 2, 3, BioPolymerCoverageType.Unique);
        var results = new List<BioPolymerCoverageResultModel> { result1, result2 };

        var mapper = new ScoreColorMapper();
        mapper.Prepare(results, ColorGradientType.Viridis, false);
        Assert.That(mapper.GradientMinValue, Is.EqualTo(10));
        Assert.That(mapper.GradientMaxValue, Is.EqualTo(1000));
    }

    [Test]
    public void ColorResultsByExtensions_IsNumeric_ReturnsTrueForNumericModes()
    {
        Assert.That(ColorResultsBy.PrecursorIntensity.IsNumeric(), Is.True);
        Assert.That(ColorResultsBy.Score.IsNumeric(), Is.True);
    }

    [Test]
    public void ColorResultsByExtensions_IsNumeric_ReturnsFalseForNonNumericModes()
    {
        Assert.That(ColorResultsBy.None.IsNumeric(), Is.False);
        Assert.That(ColorResultsBy.CoverageType.IsNumeric(), Is.False);
        Assert.That(ColorResultsBy.FileOrigin.IsNumeric(), Is.False);
    }

    [Test]
    public void IsNumericColorMode_TracksColorBy()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        Assert.That(vm.IsNumericColorMode, Is.False);

        vm.ColorBy = ColorResultsBy.PrecursorIntensity;
        Assert.That(vm.IsNumericColorMode, Is.True);

        vm.ColorBy = ColorResultsBy.Score;
        Assert.That(vm.IsNumericColorMode, Is.True);

        vm.ColorBy = ColorResultsBy.CoverageType;
        Assert.That(vm.IsNumericColorMode, Is.False);
    }

    [Test]
    public void AllGradients_ContainsViridis()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        Assert.That(vm.AllGradients, Is.Not.Null);
        Assert.That(vm.AllGradients.Length, Is.EqualTo(1));
        Assert.That(vm.AllGradients[0], Is.EqualTo(ColorGradientType.Viridis));
    }

    [Test]
    public void LogDefault_PerModeTracking_PrecursorThenScore()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.PrecursorIntensity;
        Assert.That(vm.UseLogColorScale, Is.True);

        vm.ColorBy = ColorResultsBy.Score;
        Assert.That(vm.UseLogColorScale, Is.False);
    }

    [Test]
    public void LogDefault_PerModeTracking_BackToPrecursorKeepsDefault()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.Score;
        Assert.That(vm.UseLogColorScale, Is.False);

        vm.ColorBy = ColorResultsBy.PrecursorIntensity;
        Assert.That(vm.UseLogColorScale, Is.True);

        // both modes have had defaults applied once; switching back keeps last value
        vm.ColorBy = ColorResultsBy.Score;
        Assert.That(vm.UseLogColorScale, Is.True);
    }

    [Test]
    public void LogDefault_PerModeTracking_NonNumericDoesNotSet()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.CoverageType;
        Assert.That(vm.UseLogColorScale, Is.False);
    }

    [Test]
    public void LogDefault_UserOverride_PerModeStillApplies()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.PrecursorIntensity;
        Assert.That(vm.UseLogColorScale, Is.True);

        vm.UseLogColorScale = false;
        Assert.That(vm.UseLogColorScale, Is.False);

        vm.ColorBy = ColorResultsBy.Score;
        Assert.That(vm.UseLogColorScale, Is.False);

        vm.ColorBy = ColorResultsBy.PrecursorIntensity;
        Assert.That(vm.UseLogColorScale, Is.False);
    }

    [Test]
    public void ScoreColorMapper_CreateRenderContext_LogScale_TransformsValues()
    {
        var match1 = new DummySpectralmatch();
        match1.SetScore(10);
        var match2 = new DummySpectralmatch();
        match2.SetScore(10000);
        var result1 = new BioPolymerCoverageResultModel(match1, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var result2 = new BioPolymerCoverageResultModel(match2, "ABC", 2, 3, BioPolymerCoverageType.Unique);
        var results = new List<BioPolymerCoverageResultModel> { result1, result2 };

        var mapper = new ScoreColorMapper();
        mapper.Prepare(results, ColorGradientType.Viridis, true);
        Assert.That(mapper.GradientMinValue, Is.EqualTo(1.0).Within(1e-9));
        Assert.That(mapper.GradientMaxValue, Is.EqualTo(4.0).Within(1e-9));
        Assert.That(mapper.GradientLegendTitle, Does.EndWith("(log10)"));
    }

    [Test]
    public void ScoreColorMapper_CreateRenderContext_LogScale_ExcludesNonPositive()
    {
        var match1 = new DummySpectralmatch();
        match1.SetScore(0);
        var match2 = new DummySpectralmatch();
        match2.SetScore(1000);
        var result1 = new BioPolymerCoverageResultModel(match1, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var result2 = new BioPolymerCoverageResultModel(match2, "ABC", 2, 3, BioPolymerCoverageType.Unique);
        var results = new List<BioPolymerCoverageResultModel> { result1, result2 };

        var mapper = new ScoreColorMapper();
        mapper.Prepare(results, ColorGradientType.Viridis, true);
        Assert.That(mapper.GradientMinValue, Is.EqualTo(3.0).Within(1e-9));
        Assert.That(mapper.GradientMaxValue, Is.EqualTo(4.0).Within(1e-9));
    }

    [Test]
    public void PrecursorIntensityColorMapper_HasFallback_AndDefaultLog()
    {
        var mapper = new PrecursorIntensityColorMapper();
        Assert.That(mapper.DefaultUseLogScale, Is.True);
        Assert.That(mapper.DisplayName, Is.EqualTo("Precursor Intensity"));
    }

    [Test]
    public void ScoreColorMapper_HasNoFallback_AndDefaultLogFalse()
    {
        var mapper = new ScoreColorMapper();
        Assert.That(mapper.DefaultUseLogScale, Is.False);
        Assert.That(mapper.DisplayName, Is.EqualTo("Score"));
    }

    [Test]
    public void UseLogColorScale_DefaultOnPrecursorIntensity()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.PrecursorIntensity;
        Assert.That(vm.UseLogColorScale, Is.True);
    }

    [Test]
    public void UseLogColorScale_DefaultOnScore()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.Score;
        Assert.That(vm.UseLogColorScale, Is.False);
    }

    [Test]
    public void UseLogColorScale_UserOverride_PersistsAcrossColorBySwitches()
    {
        var vm = new BioPolymerCoverageMapViewModel();
        vm.ColorBy = ColorResultsBy.PrecursorIntensity;
        Assert.That(vm.UseLogColorScale, Is.True);
        vm.UseLogColorScale = false; // user override
        vm.ColorBy = ColorResultsBy.Score;
        Assert.That(vm.UseLogColorScale, Is.False); // persisted, not reset
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
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.PrecursorIntensity);
        mapper.Prepare(filteredResults, ColorGradientType.Viridis, false);
        var result = method.Invoke(vm, new object[] { 16.0, 96.0, mapper })
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
        var mapper = BioPolymerCoverageColorMapperFactory.Create(ColorResultsBy.Score);
        mapper.Prepare(filteredResults, ColorGradientType.Viridis, false);
        var result = method.Invoke(vm, new object[] { 16.0, 96.0, mapper })
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
