using GuiFunctions.MetaDraw.BioPolymerCoverage.ColorMapping.Gradient;
using NUnit.Framework;
using System;
using System.Linq;
using System.Windows.Media;

namespace Test.MetaDraw;

[TestFixture]
public class ColorGradientTests
{
    #region Viridis Tests
    [Test]
    public void Viridis_GradientType_IsViridis()
    {
        var gradient = new ViridisColorGradient();
        Assert.That(gradient.GradientType, Is.EqualTo(ColorGradientType.Viridis));
    }

    [Test]
    public void Viridis_BinCount_Is20()
    {
        var gradient = new ViridisColorGradient();
        Assert.That(gradient.BinCount, Is.EqualTo(20));
    }

    [Test]
    public void Viridis_GetBrushes_Returns20NonNullBrushes()
    {
        var gradient = new ViridisColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes.Count, Is.EqualTo(20));
        Assert.That(brushes.All(b => b != null), Is.True);
    }

    [Test]
    public void Viridis_GetBrushes_AreUniqueInstances()
    {
        var gradient = new ViridisColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes.Distinct().Count(), Is.EqualTo(20));
    }

    [Test]
    public void Viridis_StartAndEndBrushesDiffer()
    {
        var gradient = new ViridisColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes[0].Color, Is.Not.EqualTo(brushes[brushes.Count - 1].Color));
    }

    [Test]
    public void Viridis_GetBrush_Zero_ReturnsFirstColor()
    {
        var gradient = new ViridisColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(gradient.GetBrush(0.0), Is.EqualTo(brushes[0]));
    }

    [Test]
    public void Viridis_GetBrush_One_ReturnsLastColor()
    {
        var gradient = new ViridisColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(gradient.GetBrush(1.0), Is.EqualTo(brushes[brushes.Count - 1]));
    }

    [Test]
    public void Viridis_GetBrush_Negative_ClampsToFirstColor()
    {
        var gradient = new ViridisColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(gradient.GetBrush(-0.5), Is.EqualTo(brushes[0]));
    }

    [Test]
    public void Viridis_GetBrush_AboveOne_ClampsToLastColor()
    {
        var gradient = new ViridisColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(gradient.GetBrush(1.5), Is.EqualTo(brushes[brushes.Count - 1]));
    }

    [Test]
    public void Viridis_GetBrush_MidValue_ReturnsMidBin()
    {
        var gradient = new ViridisColorGradient();
        var brushes = gradient.GetBrushes();
        // 0.5 * 19 = 9.5, clamped to 9
        Assert.That(gradient.GetBrush(0.5), Is.EqualTo(brushes[9]));
    }
    #endregion

    #region Plasma Tests
    [Test]
    public void Plasma_GradientType_IsPlasma()
    {
        var gradient = new PlasmaColorGradient();
        Assert.That(gradient.GradientType, Is.EqualTo(ColorGradientType.Plasma));
    }

    [Test]
    public void Plasma_BinCount_Is20()
    {
        var gradient = new PlasmaColorGradient();
        Assert.That(gradient.BinCount, Is.EqualTo(20));
    }

    [Test]
    public void Plasma_GetBrushes_Returns20NonNullBrushes()
    {
        var gradient = new PlasmaColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes.Count, Is.EqualTo(20));
        Assert.That(brushes.All(b => b != null), Is.True);
    }

    [Test]
    public void Plasma_StartAndEndBrushesDiffer()
    {
        var gradient = new PlasmaColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes[0].Color, Is.Not.EqualTo(brushes[brushes.Count - 1].Color));
    }
    #endregion

    #region Inferno Tests
    [Test]
    public void Inferno_GradientType_IsInferno()
    {
        var gradient = new InfernoColorGradient();
        Assert.That(gradient.GradientType, Is.EqualTo(ColorGradientType.Inferno));
    }

    [Test]
    public void Inferno_BinCount_Is20()
    {
        var gradient = new InfernoColorGradient();
        Assert.That(gradient.BinCount, Is.EqualTo(20));
    }

    [Test]
    public void Inferno_GetBrushes_Returns20NonNullBrushes()
    {
        var gradient = new InfernoColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes.Count, Is.EqualTo(20));
        Assert.That(brushes.All(b => b != null), Is.True);
    }

    [Test]
    public void Inferno_StartAndEndBrushesDiffer()
    {
        var gradient = new InfernoColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes[0].Color, Is.Not.EqualTo(brushes[brushes.Count - 1].Color));
    }
    #endregion

    #region Turbo Tests
    [Test]
    public void Turbo_GradientType_IsTurbo()
    {
        var gradient = new TurboColorGradient();
        Assert.That(gradient.GradientType, Is.EqualTo(ColorGradientType.Turbo));
    }

    [Test]
    public void Turbo_BinCount_Is20()
    {
        var gradient = new TurboColorGradient();
        Assert.That(gradient.BinCount, Is.EqualTo(20));
    }

    [Test]
    public void Turbo_GetBrushes_Returns20NonNullBrushes()
    {
        var gradient = new TurboColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes.Count, Is.EqualTo(20));
        Assert.That(brushes.All(b => b != null), Is.True);
    }

    [Test]
    public void Turbo_StartAndEndBrushesDiffer()
    {
        var gradient = new TurboColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes[0].Color, Is.Not.EqualTo(brushes[brushes.Count - 1].Color));
    }
    #endregion

    #region Grayscale Tests
    [Test]
    public void Grayscale_GradientType_IsGrayscale()
    {
        var gradient = new GrayscaleColorGradient();
        Assert.That(gradient.GradientType, Is.EqualTo(ColorGradientType.Grayscale));
    }

    [Test]
    public void Grayscale_BinCount_Is20()
    {
        var gradient = new GrayscaleColorGradient();
        Assert.That(gradient.BinCount, Is.EqualTo(20));
    }

    [Test]
    public void Grayscale_GetBrushes_Returns20NonNullBrushes()
    {
        var gradient = new GrayscaleColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes.Count, Is.EqualTo(20));
        Assert.That(brushes.All(b => b != null), Is.True);
    }

    [Test]
    public void Grayscale_StartAndEndBrushesDiffer()
    {
        var gradient = new GrayscaleColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes[0].Color, Is.Not.EqualTo(brushes[brushes.Count - 1].Color));
    }

    [Test]
    public void Grayscale_FirstBrushIsWhite()
    {
        var gradient = new GrayscaleColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes[0].Color, Is.EqualTo(Colors.White));
    }

    [Test]
    public void Grayscale_LastBrushIsBlack()
    {
        var gradient = new GrayscaleColorGradient();
        var brushes = gradient.GetBrushes();
        Assert.That(brushes[brushes.Count - 1].Color, Is.EqualTo(Colors.Black));
    }
    #endregion

    #region Parameterized Tests for All Gradients
    [TestCase(ColorGradientType.Viridis, typeof(ViridisColorGradient))]
    [TestCase(ColorGradientType.Plasma, typeof(PlasmaColorGradient))]
    [TestCase(ColorGradientType.Inferno, typeof(InfernoColorGradient))]
    [TestCase(ColorGradientType.Turbo, typeof(TurboColorGradient))]
    [TestCase(ColorGradientType.Grayscale, typeof(GrayscaleColorGradient))]
    public void AllGradients_GetBrush_Zero_ReturnsFirstColor(ColorGradientType type, Type gradientType)
    {
        var gradient = ColorGradientFactory.Create(type);
        var brushes = gradient.GetBrushes();
        Assert.That(gradient.GetBrush(0.0), Is.EqualTo(brushes[0]));
    }

    [TestCase(ColorGradientType.Viridis, typeof(ViridisColorGradient))]
    [TestCase(ColorGradientType.Plasma, typeof(PlasmaColorGradient))]
    [TestCase(ColorGradientType.Inferno, typeof(InfernoColorGradient))]
    [TestCase(ColorGradientType.Turbo, typeof(TurboColorGradient))]
    [TestCase(ColorGradientType.Grayscale, typeof(GrayscaleColorGradient))]
    public void AllGradients_GetBrush_One_ReturnsLastColor(ColorGradientType type, Type gradientType)
    {
        var gradient = ColorGradientFactory.Create(type);
        var brushes = gradient.GetBrushes();
        Assert.That(gradient.GetBrush(1.0), Is.EqualTo(brushes[brushes.Count - 1]));
    }

    [TestCase(ColorGradientType.Viridis, typeof(ViridisColorGradient))]
    [TestCase(ColorGradientType.Plasma, typeof(PlasmaColorGradient))]
    [TestCase(ColorGradientType.Inferno, typeof(InfernoColorGradient))]
    [TestCase(ColorGradientType.Turbo, typeof(TurboColorGradient))]
    [TestCase(ColorGradientType.Grayscale, typeof(GrayscaleColorGradient))]
    public void AllGradients_GetBrush_Negative_ClampsToFirstColor(ColorGradientType type, Type gradientType)
    {
        var gradient = ColorGradientFactory.Create(type);
        var brushes = gradient.GetBrushes();
        Assert.That(gradient.GetBrush(-0.5), Is.EqualTo(brushes[0]));
    }

    [TestCase(ColorGradientType.Viridis, typeof(ViridisColorGradient))]
    [TestCase(ColorGradientType.Plasma, typeof(PlasmaColorGradient))]
    [TestCase(ColorGradientType.Inferno, typeof(InfernoColorGradient))]
    [TestCase(ColorGradientType.Turbo, typeof(TurboColorGradient))]
    [TestCase(ColorGradientType.Grayscale, typeof(GrayscaleColorGradient))]
    public void AllGradients_GetBrush_AboveOne_ClampsToLastColor(ColorGradientType type, Type gradientType)
    {
        var gradient = ColorGradientFactory.Create(type);
        var brushes = gradient.GetBrushes();
        Assert.That(gradient.GetBrush(1.5), Is.EqualTo(brushes[brushes.Count - 1]));
    }

    [TestCase(ColorGradientType.Viridis, typeof(ViridisColorGradient))]
    [TestCase(ColorGradientType.Plasma, typeof(PlasmaColorGradient))]
    [TestCase(ColorGradientType.Inferno, typeof(InfernoColorGradient))]
    [TestCase(ColorGradientType.Turbo, typeof(TurboColorGradient))]
    [TestCase(ColorGradientType.Grayscale, typeof(GrayscaleColorGradient))]
    public void AllGradients_GetBrush_MidValue_ReturnsMidBin(ColorGradientType type, Type gradientType)
    {
        var gradient = ColorGradientFactory.Create(type);
        var brushes = gradient.GetBrushes();
        // 0.5 * 19 = 9.5, clamped to 9
        Assert.That(gradient.GetBrush(0.5), Is.EqualTo(brushes[9]));
    }
    #endregion

    #region Factory Tests
    [Test]
    public void ColorGradientFactory_Create_Viridis_ReturnsViridisColorGradient()
    {
        var gradient = ColorGradientFactory.Create(ColorGradientType.Viridis);
        Assert.That(gradient, Is.InstanceOf<ViridisColorGradient>());
    }

    [Test]
    public void ColorGradientFactory_Create_Plasma_ReturnsPlasmaColorGradient()
    {
        var gradient = ColorGradientFactory.Create(ColorGradientType.Plasma);
        Assert.That(gradient, Is.InstanceOf<PlasmaColorGradient>());
    }

    [Test]
    public void ColorGradientFactory_Create_Inferno_ReturnsInfernoColorGradient()
    {
        var gradient = ColorGradientFactory.Create(ColorGradientType.Inferno);
        Assert.That(gradient, Is.InstanceOf<InfernoColorGradient>());
    }

    [Test]
    public void ColorGradientFactory_Create_Turbo_ReturnsTurboColorGradient()
    {
        var gradient = ColorGradientFactory.Create(ColorGradientType.Turbo);
        Assert.That(gradient, Is.InstanceOf<TurboColorGradient>());
    }

    [Test]
    public void ColorGradientFactory_Create_Grayscale_ReturnsGrayscaleColorGradient()
    {
        var gradient = ColorGradientFactory.Create(ColorGradientType.Grayscale);
        Assert.That(gradient, Is.InstanceOf<GrayscaleColorGradient>());
    }

    [Test]
    public void ColorGradientFactory_Create_Unknown_Throws()
    {
        var unknown = (ColorGradientType)9999;
        Assert.Throws<ArgumentOutOfRangeException>(() => ColorGradientFactory.Create(unknown));
    }
    #endregion
}
