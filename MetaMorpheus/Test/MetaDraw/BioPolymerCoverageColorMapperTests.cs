using GuiFunctions.MetaDraw;
using NUnit.Framework;
using System;
using System.Windows.Media;

namespace Test.MetaDraw;

[TestFixture]
public class BioPolymerCoverageColorMapperTests
{
    private static BioPolymerCoverageResultModel MakeResult(double? precursorIntensity, double score)
    {
        var psm = new DummySpectralmatch();
        psm.SetScore(score);
        if (precursorIntensity.HasValue)
        {
            var backing = typeof(Readers.SpectrumMatchFromTsv)
                .GetField("<PrecursorIntensity>k__BackingField",
                    System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
            backing.SetValue(psm, precursorIntensity.Value);
        }
        return new BioPolymerCoverageResultModel(psm, "ABC", 1, 2, BioPolymerCoverageType.Unique);
    }

    #region Categorical

    [Test]
    public void CategoricalMapper_IsNotNumeric()
    {
        var mapper = new CategoricalBioPolymerCoverageColorMapper(
            ColorResultsBy.CoverageType,
            r => Brushes.Black);
        Assert.That(mapper.IsNumeric, Is.False);
    }

    [Test]
    public void CategoricalMapper_ColorBy_ReturnsPassedEnum()
    {
        var mapper = new CategoricalBioPolymerCoverageColorMapper(
            ColorResultsBy.CoverageType,
            r => Brushes.Black);
        Assert.That(mapper.ColorBy, Is.EqualTo(ColorResultsBy.CoverageType));
    }

    [Test]
    public void CategoricalMapper_GetBrush_UsesSelector()
    {
        var mapper = new CategoricalBioPolymerCoverageColorMapper(
            ColorResultsBy.CoverageType,
            r => Brushes.Red);
        var result = MakeResult(100, 10);
        mapper.Prepare([result], ColorGradientType.Viridis, false);
        Assert.That(mapper.GetBrush(result), Is.EqualTo(Brushes.Red));
    }

    [Test]
    public void CategoricalMapper_Constructor_ThrowsOnNullSelector()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new CategoricalBioPolymerCoverageColorMapper(ColorResultsBy.CoverageType, null!));
    }

    [Test]
    public void CategoricalMapper_GetBrush_ScaleIsIgnored()
    {
        var mapper = new CategoricalBioPolymerCoverageColorMapper(
            ColorResultsBy.CoverageType,
            r => Brushes.Blue);
        var result = MakeResult(null, 5);
        mapper.Prepare([result], ColorGradientType.Viridis, false);
        Assert.That(mapper.GetBrush(result), Is.EqualTo(Brushes.Blue));
    }

    #endregion

    #region Numeric PrecursorIntensity

    [Test]
    public void PrecursorIntensityMapper_IsNumeric()
    {
        var mapper = new PrecursorIntensityColorMapper();
        Assert.That(mapper.IsNumeric, Is.True);
    }

    [Test]
    public void PrecursorIntensityMapper_ColorBy_IsPrecursorIntensity()
    {
        var mapper = new PrecursorIntensityColorMapper();
        Assert.That(mapper.ColorBy, Is.EqualTo(ColorResultsBy.PrecursorIntensity));
    }

    [Test]
    public void PrecursorIntensityMapper_GetNumericValue_ReturnsValue()
    {
        var result = MakeResult(500000, 10);
        var mapper = new PrecursorIntensityColorMapper();
        Assert.That(mapper.GetNumericValue(result), Is.EqualTo(500000));
    }

    [Test]
    public void PrecursorIntensityMapper_GetBrush_MinValue_ReturnsFirstColor()
    {
        var mapper = new PrecursorIntensityColorMapper();
        var result = MakeResult(100, 0);
        var maxResult = MakeResult(1000, 0);
        mapper.Prepare([result, maxResult], ColorGradientType.Viridis, false);
        Assert.That(mapper.GetBrush(result), Is.EqualTo(new ViridisColorGradient().GetBrushes()[0]));
    }

    [Test]
    public void PrecursorIntensityMapper_GetBrush_MaxValue_ReturnsLastColor()
    {
        var mapper = new PrecursorIntensityColorMapper();
        var minResult = MakeResult(100, 0);
        var result = MakeResult(1000, 0);
        mapper.Prepare([minResult, result], ColorGradientType.Viridis, false);
        Assert.That(mapper.GetBrush(result), Is.EqualTo(new ViridisColorGradient().GetBrushes()[19]));
    }

    [Test]
    public void PrecursorIntensityMapper_GetBrush_NullValue_ReturnsGray()
    {
        var psm = new DummySpectralmatch();
        psm.SetScore(50);
        var backing = typeof(Readers.SpectrumMatchFromTsv)
            .GetField("<PrecursorIntensity>k__BackingField",
                System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance);
        backing.SetValue(psm, null);
        var result = new BioPolymerCoverageResultModel(psm, "ABC", 1, 2, BioPolymerCoverageType.Unique);
        var mapper = new PrecursorIntensityColorMapper();
        mapper.Prepare([result], ColorGradientType.Viridis, false);
        var brush = mapper.GetBrush(result);
        Assert.That(brush.Color, Is.EqualTo(Colors.Gray));
    }

    [Test]
    public void PrecursorIntensityMapper_GetBrush_EmptyPreparedData_UsesFirstGradientColor()
    {
        var mapper = new PrecursorIntensityColorMapper();
        var result = MakeResult(500, 0);
        mapper.Prepare([], ColorGradientType.Viridis, false);
        var brush = mapper.GetBrush(result);
        Assert.That(brush, Is.EqualTo(new ViridisColorGradient().GetBrushes()[0]));
    }

    #endregion

    #region Numeric Score

    [Test]
    public void ScoreMapper_IsNumeric()
    {
        var mapper = new ScoreColorMapper();
        Assert.That(mapper.IsNumeric, Is.True);
    }

    [Test]
    public void ScoreMapper_ColorBy_IsScore()
    {
        var mapper = new ScoreColorMapper();
        Assert.That(mapper.ColorBy, Is.EqualTo(ColorResultsBy.Score));
    }

    [Test]
    public void ScoreMapper_GetNumericValue_ReturnsScore()
    {
        var result = MakeResult(500, 42.5);
        var mapper = new ScoreColorMapper();
        Assert.That(mapper.GetNumericValue(result), Is.EqualTo(42.5));
    }

    [Test]
    public void ScoreMapper_GetBrush_MidValue_UsesNormalizedGradient()
    {
        var mapper = new ScoreColorMapper();
        var result = MakeResult(null, 50);
        mapper.Prepare([MakeResult(null, 0), result, MakeResult(null, 100)], ColorGradientType.Viridis, false);
        var expectedBin = (int)Math.Clamp(0.5 * 19, 0, 19);
        Assert.That(mapper.GetBrush(result), Is.EqualTo(new ViridisColorGradient().GetBrushes()[expectedBin]));
    }

    [Test]
    public void ScoreMapper_GetBrush_NullResult_ReturnsGray()
    {
        var mapper = new ScoreColorMapper();
        mapper.Prepare([], ColorGradientType.Viridis, false);
        var brush = mapper.GetBrush(null!);
        Assert.That(brush.Color, Is.EqualTo(Colors.Gray));
    }

    #endregion

    #region Legend title

    [Test]
    public void PrecursorIntensityMapper_LegendTitle_IsDisplayName()
    {
        var mapper = new PrecursorIntensityColorMapper();
        mapper.Prepare([], ColorGradientType.Viridis, false);
        Assert.That(mapper.GradientLegendTitle, Is.EqualTo("Precursor Intensity"));
    }

    [Test]
    public void CategoricalMapper_LegendTitle_IsEnumName()
    {
        var mapper = new CategoricalBioPolymerCoverageColorMapper(
            ColorResultsBy.CoverageType,
            r => Brushes.Black);
        mapper.Prepare([], ColorGradientType.Viridis, false);
        Assert.That(mapper.LegendItems, Is.Empty);
    }

    #endregion
}
