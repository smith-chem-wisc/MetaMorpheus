using GuiFunctions.MetaDraw;
using NUnit.Framework;
using System;
using System.Linq;
using System.Windows.Media;

namespace Test.MetaDraw;

[TestFixture]
public class ColorGradientTests
{
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

    [Test]
    public void ColorGradientFactory_Create_Viridis_ReturnsViridisColorGradient()
    {
        var gradient = ColorGradientFactory.Create(ColorGradientType.Viridis);
        Assert.That(gradient, Is.InstanceOf<ViridisColorGradient>());
    }

    [Test]
    public void ColorGradientFactory_Create_Unknown_Throws()
    {
        var unknown = (ColorGradientType)9999;
        Assert.Throws<ArgumentOutOfRangeException>(() => ColorGradientFactory.Create(unknown));
    }
}
