using System.Collections.Generic;
using System.Linq;
using EngineLayer.FdrAnalysis;
using EngineLayer.ParallelSearch.FdrAlignment;
using NUnit.Framework;

namespace Test.ParallelSearchTask.FdrAlignment;

/// <summary>
/// A transient match takes the FDR values of the baseline entry at its score. The PEP q-value is not monotone in score,
/// so it used to be whatever one neighbouring baseline PSM had; the lookups now hold its isotonic fit.
/// </summary>
[TestFixture]
public class MonotonePepQValueTests
{
    private static (double Score, FdrInfo Value) Entry(double score, double pepQ) =>
        (score, new FdrInfo { PEP_QValue = pepQ });

    [Test]
    public void TheFit_NeverFalls_AsTheScoreFalls_AndAveragesViolators()
    {
        var entries = new List<(double Score, FdrInfo Value)>
        {
            Entry(30, 0.001),
            Entry(25, 0.20), // a neighbour with a bad PEP
            Entry(24, 0.002),
            Entry(20, 0.01),
            Entry(15, 0.30),
        };

        var fitted = BaselineFdr.WithMonotonePepQValues(entries)
            .OrderByDescending(e => e.Score)
            .Select(e => e.Value.PEP_QValue)
            .ToList();

        for (int i = 1; i < fitted.Count; i++)
            Assert.That(fitted[i], Is.GreaterThanOrEqualTo(fitted[i - 1]), "non-decreasing as the score falls");
        // 0.20, 0.002 and 0.01 pool to their mean
        Assert.That(fitted[1], Is.EqualTo((0.20 + 0.002 + 0.01) / 3).Within(1e-12));
        Assert.That(fitted[0], Is.EqualTo(0.001));
        Assert.That(fitted[4], Is.EqualTo(0.30));
    }

    [Test]
    public void UnsetPepQValues_AreLeftOut_AndLeftAlone()
    {
        var entries = new List<(double Score, FdrInfo Value)> { Entry(30, 0.05), Entry(25, 2), Entry(20, 0.01) };

        var fitted = BaselineFdr.WithMonotonePepQValues(entries).ToDictionary(e => e.Score, e => e.Value.PEP_QValue);

        Assert.That(fitted[25], Is.EqualTo(2), "unset stays unset");
        Assert.That(fitted[30], Is.EqualTo(0.03).Within(1e-12), "0.05 above 0.01 pools to 0.03");
        Assert.That(fitted[20], Is.EqualTo(0.03).Within(1e-12));
    }

    [Test]
    public void EqualScores_ShareOneValue()
    {
        var entries = new List<(double Score, FdrInfo Value)> { Entry(20, 0.01), Entry(20, 0.03), Entry(10, 0.5) };

        var fitted = BaselineFdr.WithMonotonePepQValues(entries);

        Assert.That(fitted.Where(e => e.Score == 20).Select(e => e.Value.PEP_QValue), Is.All.EqualTo(0.02).Within(1e-12));
    }
}
