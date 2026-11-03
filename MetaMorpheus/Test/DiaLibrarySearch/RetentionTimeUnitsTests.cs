#nullable enable
using System;
using System.Diagnostics.CodeAnalysis;
using EngineLayer.DiaLibrarySearch;
using NUnit.Framework;

namespace Test.DiaLibrarySearch;

[TestFixture]
[ExcludeFromCodeCoverage]
public class RetentionTimeUnitsTests
{
    [Test]
    public void EachRetentionTimeNamesItsUnitWhenPrinted()
    {
        Assert.That(new Irt(12).ToString(), Is.EqualTo("12 iRT"));
        Assert.That(new RtMinutes(3).ToString(), Is.EqualTo("3 min"));
    }

    [Test]
    public void ALinearMapInvertsExactly()
    {
        var map = LinearIrtMap.ThroughPoints((new RtMinutes(1), new Irt(-20)), (new RtMinutes(9), new Irt(120)));

        Assert.That(map.ToIrt(new RtMinutes(1)).Value, Is.EqualTo(-20).Within(1e-9));
        Assert.That(map.ToRtMinutes(new Irt(120)).Value, Is.EqualTo(9).Within(1e-9));
        foreach (double minutes in new[] { 0, 1, 4.2, 9, 12 })
            Assert.That(map.ToRtMinutes(map.ToIrt(new RtMinutes(minutes))).Value, Is.EqualTo(minutes).Within(1e-9), $"{minutes} min");
    }

    /// <summary>A slope of zero, or one that is not finite, cannot be inverted, so the map would have no way back.</summary>
    [TestCase(0)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    public void ASlopeThatCannotBeInvertedIsRefused(double slope)
    {
        var e = Assert.Throws<ArgumentOutOfRangeException>(() => new LinearIrtMap(slope, 0));
        Assert.That(e!.ParamName, Is.EqualTo("slope"));
    }

    [TestCase(double.NaN)]
    [TestCase(double.NegativeInfinity)]
    public void AnInterceptThatIsNotFiniteIsRefused(double intercept)
    {
        var e = Assert.Throws<ArgumentOutOfRangeException>(() => new LinearIrtMap(1, intercept));
        Assert.That(e!.ParamName, Is.EqualTo("intercept"));
    }

    /// <summary>Two anchors at one retention time give an infinite slope, and no map.</summary>
    [Test]
    public void TwoAnchorsAtOneRetentionTimeAreRefused()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            LinearIrtMap.ThroughPoints((new RtMinutes(4), new Irt(0)), (new RtMinutes(4), new Irt(50))));
    }
}
