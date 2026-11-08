#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer.DiaLibrarySearch;
using NUnit.Framework;

namespace Test.DiaLibrarySearch;

/// <summary>
/// A run's MS1 mass offset as a function of retention time, fitted from confident identifications' signed MS1 errors: the
/// median error in equal-count RT bins, interpolated between bin centres and held flat beyond them. On PXD022589 the offset
/// drifts from +0.1 to +2.9 ppm over the gradient, so a fixed window centred at 0 clips real precursors unevenly.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class Ms1OffsetModelTests
{
    [Test]
    public void AnOffsetThatDriftsWithRetentionTimeIsRecovered()
    {
        var random = new Random(3);
        var points = Enumerable.Range(0, 2000).Select(i =>
        {
            double rt = 60.0 * i / 2000;
            return (rt, 1 + 0.1 * rt + (random.NextDouble() - 0.5));
        }).ToList();

        var model = Ms1OffsetModel.Fit(points);

        foreach (double rt in new[] { 5.0, 20.0, 40.0, 55.0 })
            Assert.That(model.OffsetPpm(rt), Is.EqualTo(1 + 0.1 * rt).Within(0.4), $"at {rt} min");
        Assert.That(model.OffsetPpm(-10), Is.EqualTo(model.OffsetPpm(0)).Within(1e-9), "flat before the first bin");
        Assert.That(model.OffsetPpm(500), Is.EqualTo(model.OffsetPpm(60)).Within(1e-9), "flat after the last bin");
    }

    [Test]
    public void FewOrNoPointsGiveAConstantOrZeroOffset()
    {
        Assert.That(Ms1OffsetModel.Fit([]).OffsetPpm(10), Is.EqualTo(0));
        var few = new List<(double, double)> { (1, 2.0), (2, 3.0), (3, 2.5) };
        Assert.That(Ms1OffsetModel.Fit(few).OffsetPpm(50), Is.EqualTo(2.5).Within(1e-9), "one bin: the median");
        Assert.That(Ms1OffsetModel.Constant(3).OffsetPpm(12), Is.EqualTo(3));
        Assert.That(Ms1OffsetModel.Fit([(1.0, double.NaN), (2.0, 4.0)]).OffsetPpm(1), Is.EqualTo(4), "non-finite errors are ignored");
    }

    /// <summary>
    /// The spread of the errors around the fitted offset, robustly (1.4826 x the median absolute residual, an SD for normal
    /// errors): what is left once the drift is removed, so a run's MS1 tolerance can follow its own mass accuracy. Errors
    /// alternating 1 ppm either side of a drifting offset give 1.4826; a constant model has no spread (NaN).
    /// </summary>
    [Test]
    public void TheSpreadAroundTheOffsetIsMeasured()
    {
        var points = Enumerable.Range(0, 2000).Select(i =>
        {
            double rt = 60.0 * i / 2000;
            return (rt, 1 + 0.1 * rt + (i % 2 == 0 ? 1.0 : -1.0));
        }).ToList();

        Assert.That(Ms1OffsetModel.Fit(points).ResidualSpreadPpm, Is.EqualTo(1.4826).Within(0.05));
        Assert.That(Ms1OffsetModel.Constant(3).ResidualSpreadPpm, Is.NaN);
        Assert.That(Ms1OffsetModel.Fit([]).ResidualSpreadPpm, Is.NaN);
    }
}
