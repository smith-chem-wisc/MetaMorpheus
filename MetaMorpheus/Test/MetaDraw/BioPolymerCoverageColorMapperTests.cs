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
        Assert.That(mapper.GetBrush(result, null), Is.EqualTo(Brushes.Red));
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
        var gradient = new ViridisColorGradient();
        var scale = new BioPolymerCoverageColorScale(0, 100, gradient);
        var mapper = new CategoricalBioPolymerCoverageColorMapper(
            ColorResultsBy.CoverageType,
            r => Brushes.Blue);
        var result = MakeResult(null, 5);
        Assert.That(mapper.GetBrush(result, scale), Is.EqualTo(Brushes.Blue));
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
        var gradient = new ViridisColorGradient();
        var scale = new BioPolymerCoverageColorScale(100, 1000, gradient);
        var mapper = new PrecursorIntensityColorMapper();
        var result = MakeResult(100, 0);
        Assert.That(mapper.GetBrush(result, scale), Is.EqualTo(gradient.GetBrushes()[0]));
    }

    [Test]
    public void PrecursorIntensityMapper_GetBrush_MaxValue_ReturnsLastColor()
    {
        var gradient = new ViridisColorGradient();
        var scale = new BioPolymerCoverageColorScale(100, 1000, gradient);
        var mapper = new PrecursorIntensityColorMapper();
        var result = MakeResult(1000, 0);
        Assert.That(mapper.GetBrush(result, scale), Is.EqualTo(gradient.GetBrushes()[19]));
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
        var brush = mapper.GetBrush(result, new BioPolymerCoverageColorScale(0, 100, new ViridisColorGradient()));
        Assert.That(brush.Color, Is.EqualTo(Colors.Gray));
    }

    [Test]
    public void PrecursorIntensityMapper_GetBrush_NullScale_ReturnsGray()
    {
        var mapper = new PrecursorIntensityColorMapper();
        var result = MakeResult(500, 0);
        var brush = mapper.GetBrush(result, null);
        Assert.That(brush.Color, Is.EqualTo(Colors.Gray));
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
        var gradient = new ViridisColorGradient();
        var scale = new BioPolymerCoverageColorScale(0, 100, gradient);
        var mapper = new ScoreColorMapper();
        var result = MakeResult(null, 50);
        var expectedBin = (int)Math.Clamp(0.5 * 19, 0, 19);
        Assert.That(mapper.GetBrush(result, scale), Is.EqualTo(gradient.GetBrushes()[expectedBin]));
    }

    [Test]
    public void ScoreMapper_GetBrush_NullResult_ReturnsGray()
    {
        var mapper = new ScoreColorMapper();
        var brush = mapper.GetBrush(null!, new BioPolymerCoverageColorScale(0, 100, new ViridisColorGradient()));
        Assert.That(brush.Color, Is.EqualTo(Colors.Gray));
    }

    #endregion

    #region Legend title

    [Test]
    public void PrecursorIntensityMapper_LegendTitle_IsEnumName()
    {
        var mapper = new PrecursorIntensityColorMapper();
        Assert.That(mapper.GetLegendTitle(null), Is.EqualTo(nameof(ColorResultsBy.PrecursorIntensity)));
    }

    [Test]
    public void CategoricalMapper_LegendTitle_IsEnumName()
    {
        var mapper = new CategoricalBioPolymerCoverageColorMapper(
            ColorResultsBy.CoverageType,
            r => Brushes.Black);
        Assert.That(mapper.GetLegendTitle(null), Is.EqualTo(nameof(ColorResultsBy.CoverageType)));
    }

    #endregion
}
