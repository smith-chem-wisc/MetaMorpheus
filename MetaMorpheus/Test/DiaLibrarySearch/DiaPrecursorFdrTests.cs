#nullable enable
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer.DiaLibrarySearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Test.DiaLibrarySearch;

/// <summary>
/// The global precursor q-value (PLAN D4): each precursor's best score across the runs, then target-decoy competition, by
/// the same q-value rule each run uses. With one run it is that run's q-value.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaPrecursorFdrTests
{
    private static DiaPrecursorMatch Match(int precursor, bool decoy, double score, double q = double.NaN) =>
        new(precursor, "PEPTIDE" + precursor, 2, 400, decoy, new Irt(0), new RtMinutes(1), new Irt(0), score, [], q);

    [Test]
    public void EachPrecursorCompetesOnceWithItsBestScore()
    {
        var run1 = new List<DiaPrecursorMatch> { Match(1, false, 0.9), Match(2, false, 0.2), Match(10, true, 0.5) };
        var run2 = new List<DiaPrecursorMatch> { Match(1, false, 0.1), Match(2, false, 0.8), Match(3, false, 0.3), Match(10, true, 0.6), Match(11, true, 0.25) };

        var global = DiaPrecursorFdr.GlobalQValues([run1, run2]);

        // Best per precursor: targets 1 (0.9), 2 (0.8), 3 (0.3); decoys 10 (0.6), 11 (0.25)
        double[] expected = DeconvolutionQValueCalculator.AssignQValues([0.9, 0.8, 0.3], [0.6, 0.25]);
        Assert.That(global.Keys, Is.EquivalentTo(new[] { 1, 2, 3 }), "targets only");
        Assert.That(global[1], Is.EqualTo(expected[0]));
        Assert.That(global[2], Is.EqualTo(expected[1]));
        Assert.That(global[3], Is.EqualTo(expected[2]));
        Assert.That(global[3], Is.GreaterThan(global[1]), "the decoy at 0.6 outranks target 3");
    }

    [Test]
    public void WithOneRunTheGlobalQValueIsTheRuns()
    {
        var targets = new[] { 0.9, 0.7, 0.4, 0.2 };
        var decoys = new[] { 0.5, 0.1 };
        double[] runQ = DeconvolutionQValueCalculator.AssignQValues(targets, decoys);
        var run = targets.Select((s, i) => Match(i, false, s, runQ[i])).Concat(decoys.Select((s, i) => Match(100 + i, true, s))).ToList();

        var global = DiaPrecursorFdr.GlobalQValues([run]);

        foreach (var m in run.Where(m => !m.IsDecoy))
            Assert.That(global[m.PrecursorIndex], Is.EqualTo(m.QValue));
    }

    /// <summary>With no decoy scored in any run, D is 0 and (D+1)/T gives every target 1/T, as each run's own q-values do; never 0.</summary>
    [Test]
    public void WithNoDecoyScoredEveryTargetGetsOneOverT()
    {
        var run1 = new List<DiaPrecursorMatch> { Match(1, false, 0.9), Match(2, false, 0.2) };
        var run2 = new List<DiaPrecursorMatch> { Match(2, false, 0.8), Match(3, false, 0.3), Match(4, false, 0.1) };

        var global = DiaPrecursorFdr.GlobalQValues([run1, run2]);

        Assert.That(global.Keys, Is.EquivalentTo(new[] { 1, 2, 3, 4 }));
        Assert.That(global.Values, Is.All.EqualTo(1.0 / 4));
    }
}
