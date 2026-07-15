using GuiFunctions.MetaDraw;
using GuiFunctions.MetaDraw.BioPolymerCoverage.ColorMapping.Gradient;
using NUnit.Framework;
using System;
using System.Linq;
using System.Windows.Media;

namespace Test.MetaDraw;

[TestFixture]
public class BioPolymerCoverageColorScaleTests
{
    [Test]
    public void Constructor_StoresMinMaxAndGradient()
    {
        var gradient = new ViridisColorGradient();
        var scale = new BioPolymerCoverageColorScale(0, 100, gradient);
        Assert.That(scale.MinValue, Is.EqualTo(0));
        Assert.That(scale.MaxValue, Is.EqualTo(100));
        Assert.That(scale.Gradient, Is.SameAs(gradient));
    }

    [Test]
    public void Constructor_Throws_OnNullGradient()
    {
        Assert.Throws<ArgumentNullException>(() => new BioPolymerCoverageColorScale(0, 100, null!));
    }

    [Test]
    public void Range_ComputesMaxMinusMin()
    {
        var scale = new BioPolymerCoverageColorScale(10, 30, new ViridisColorGradient());
        Assert.That(scale.Range, Is.EqualTo(20));
    }

    [Test]
    public void Normalize_MinValue_ReturnsZero()
    {
        var scale = new BioPolymerCoverageColorScale(0, 100, new ViridisColorGradient());
        Assert.That(scale.Normalize(0), Is.EqualTo(0));
    }

    [Test]
    public void Normalize_MaxValue_ReturnsOne()
    {
        var scale = new BioPolymerCoverageColorScale(0, 100, new ViridisColorGradient());
        Assert.That(scale.Normalize(100), Is.EqualTo(1));
    }

    [Test]
    public void Normalize_MidValue_ReturnsHalf()
    {
        var scale = new BioPolymerCoverageColorScale(0, 100, new ViridisColorGradient());
        Assert.That(scale.Normalize(50), Is.EqualTo(0.5));
    }

    [Test]
    public void Normalize_BelowMin_ClampsToZero()
    {
        var scale = new BioPolymerCoverageColorScale(0, 100, new ViridisColorGradient());
        Assert.That(scale.Normalize(-10), Is.EqualTo(0));
    }

    [Test]
    public void Normalize_AboveMax_ClampsToOne()
    {
        var scale = new BioPolymerCoverageColorScale(0, 100, new ViridisColorGradient());
        Assert.That(scale.Normalize(150), Is.EqualTo(1));
    }

    [Test]
    public void Normalize_ZeroRange_ReturnsZero()
    {
        var scale = new BioPolymerCoverageColorScale(5, 5, new ViridisColorGradient());
        Assert.That(scale.Normalize(5), Is.EqualTo(0));
    }
}
