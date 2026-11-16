#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer.DiaLibrarySearch;
using NUnit.Framework;

namespace Test.DiaLibrarySearch;

/// <summary>
/// What a run's confident targets say of the library's iRT, by sequence composition. A predicted iRT is fixed for some other column, and on both benchmark runs what is left after the run's RT->iRT
/// curve has a shape: K, R and H elute 1.5 to 3 iRT earlier than predicted. A fit on 1,000 anchors takes the median error of the targets a search scores too low from 2.47 to 1.59 and from 2.22 to 1.25.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaIrtCompositionModelTests
{
    private const string Acids = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>Random tryptic-looking peptides, each with the residual the rule gives it.</summary>
    private static List<DiaIrtCompositionModel.Anchor> Anchors(int count, Func<string, double> residual, int seed = 1)
    {
        var random = new Random(seed);
        return Enumerable.Range(0, count).Select(_ =>
        {
            int length = random.Next(7, 25);
            string sequence = new string(Enumerable.Range(0, length - 1).Select(_ => Acids[random.Next(Acids.Length)]).ToArray()) + (random.Next(2) == 0 ? 'K' : 'R');
            return new DiaIrtCompositionModel.Anchor(sequence, residual(sequence));
        }).ToList();
    }

    private static int Count(string sequence, char residue) => sequence.Count(c => c == residue);

    [Test]
    public void ItLearnsWhatEachResidueDoes()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(2000, s => -3 * Count(s, 'K') + 2 * Count(s, 'P')));

        Assert.That(model.Offset("AAAAAAAAKK") - model.Offset("AAAAAAAAAK"), Is.EqualTo(-3).Within(0.2));
        Assert.That(model.Offset("APAAAAAAAK") - model.Offset("AAAAAAAAAK"), Is.EqualTo(2).Within(0.2));
    }

    /// <summary>A permuted decoy shares its target's composition, so it moves as its target does, and the target-decoy estimate stays fair.</summary>
    [Test]
    public void APermutationOfTheSameCompositionMovesAlike()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(2000, s => -3 * Count(s, 'H') + Count(s, 'V')));

        Assert.That(model.Offset("AVGHLTEK"), Is.EqualTo(model.Offset("AGVLHTEK")).Within(1e-9));
    }

    /// <summary>Short peptides read differently from what their residues add up to.</summary>
    [Test]
    public void ItLearnsWhatLengthDoes()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(4000, s => 40.0 / s.Length));
        string shortOne = new string('A', 7) + "K", longOne = new string('A', 19) + "K";

        // Against the residues' own contribution, which is the same term by term, what is left is the length
        Assert.That(model.Offset(shortOne) - model.Offset(longOne), Is.EqualTo(40.0 / 8 - 40.0 / 20).Within(0.5));
    }

    /// <summary>Basic and acidic residues at the N-terminus do not read as they do further in: the free amine is right beside them.</summary>
    [TestCase("KRH", 'K')]
    [TestCase("DE", 'D')]
    public void ItLearnsWhatTheFirstResidueDoes(string group, char residue)
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(4000, s => group.Contains(s[0]) ? 3 : 0));
        string first = residue + "AAAAAAAAK", later = "AAAAA" + residue + "AAAK";

        Assert.That(model.Offset(first) - model.Offset(later), Is.EqualTo(3).Within(0.6));
    }

    [Test]
    public void ItLearnsWhatABasicSecondResidueDoes()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(4000, s => "KRH".Contains(s[1]) ? 3 : 0));

        Assert.That(model.Offset("ARAAAAAAAK") - model.Offset("AAAAARAAAK"), Is.EqualTo(3).Within(0.6));
    }

    /// <summary>It only says what is left after the run's curve, so a run with nothing left moves nothing.</summary>
    [Test]
    public void WithNothingLeftNothingMoves()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(1000, _ => 0));

        Assert.That(model.Offset("AVGHLTEK"), Is.EqualTo(0).Within(1e-9));
    }

    [Test]
    public void TheMeanOfWhatIsLeftIsLearnedToo()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(1000, _ => 1.5));

        Assert.That(model.Offset("AVGHLTEK"), Is.EqualTo(1.5).Within(0.1));
        Assert.That(model.Offset("AVGHLTEKAVGHLTEK"), Is.EqualTo(1.5).Within(0.1));
        Assert.That(model.Offset("XXXXXXXXXX"), Is.EqualTo(1.5).Within(0.2), "a constant is not something the residues add up to");
    }

    [Test]
    public void TooFewAnchorsMoveNothing()
    {
        int below = DiaIrtCompositionModel.MinimumAnchors - 1;

        Assert.That(DiaIrtCompositionModel.Fit(Anchors(below, s => 5 * Count(s, 'K'))).Offset("AAAAKK"), Is.EqualTo(0));
        Assert.That(DiaIrtCompositionModel.Fit(Anchors(below + 1, s => 5 * Count(s, 'K'))).Offset("AAAAKK"), Is.Not.EqualTo(0));
        Assert.That(DiaIrtCompositionModel.MinimumAnchors, Is.EqualTo(200));
    }

    /// <summary>A residue no anchor holds has nothing to say, and a library holding it must not be thrown far off. That also takes a penalty, or the fit could not be solved.</summary>
    [Test]
    public void AResidueNoAnchorHoldsDoesNothing()
    {
        var anchors = Anchors(2000, s => -2 * Count(s, 'K')).Where(a => !a.BaseSequence.Contains('W')).ToList();
        var model = DiaIrtCompositionModel.Fit(anchors);

        Assert.That(model.Offset("AAAWAAAK") - model.Offset("AAAAAAAK"), Is.EqualTo(0).Within(0.2), "it has one A fewer, and no term of its own for W");
    }

    [Test]
    public void AnUnknownResidueIsIgnored()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(2000, s => Count(s, 'K')));

        Assert.That(double.IsFinite(model.Offset("XXXXXXXK")), Is.True);
        Assert.That(model.Offset("XXXXXXXK") - model.Offset("XXXXXXXX"), Is.EqualTo(1).Within(0.2), "an unknown residue counts as nothing, and the K still does");
    }

    /// <summary>A calibration's anchors are confident, but some sit on the wrong peak. In one pass, these would pull the residue they hold well off.</summary>
    [Test]
    public void AnchorsOnTheWrongPeakDoNotMoveItFar()
    {
        var anchors = Anchors(2000, s => -2 * Count(s, 'K')).Select((a, i) => i % 20 == 0 && Count(a.BaseSequence, 'K') >= 2 ? a with { Residual = 60 } : a).ToList();
        var model = DiaIrtCompositionModel.Fit(anchors);

        Assert.That(model.Offset("AAAAAAAAKK") - model.Offset("AAAAAAAAAK"), Is.EqualTo(-2).Within(0.2));
    }

    // Calibration fitting the model: DiaIrtSelfCalibrationTests

    [Test]
    public void AnEmptySequenceDoesNotThrow()
    {
        var model = DiaIrtCompositionModel.Fit(Anchors(1000, s => Count(s, 'K')));

        Assert.That(double.IsFinite(model.Offset("")), Is.True);
        Assert.That(double.IsFinite(model.Offset("K")), Is.True);
    }
}
