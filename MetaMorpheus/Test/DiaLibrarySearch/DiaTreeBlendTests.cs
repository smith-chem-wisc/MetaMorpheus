#nullable enable
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using EngineLayer.DiaLibrarySearch;
using NUnit.Framework;

namespace Test.DiaLibrarySearch;

[TestFixture]
[ExcludeFromCodeCoverage]
public class DiaTreeBlendTests
{
    private static double Gaussian(Random random) =>
        Math.Sqrt(-2 * Math.Log(1 - random.NextDouble())) * Math.Cos(2 * Math.PI * random.NextDouble());

    /// <summary>The chance a target outranks a decoy, ties counting half.</summary>
    private static double Auc(IReadOnlyList<double> scores, IReadOnlyList<bool> isDecoy)
    {
        var targets = scores.Where((_, i) => !isDecoy[i]).ToArray();
        var decoys = scores.Where((_, i) => isDecoy[i]).Order().ToArray();
        double wins = 0;
        foreach (double target in targets)
        {
            int below = Array.BinarySearch(decoys, target);
            if (below < 0)
            {
                wins += ~below;
                continue;
            }
            int low = below, high = below;
            while (low > 0 && decoys[low - 1] == target) low--;
            while (high + 1 < decoys.Length && decoys[high + 1] == target) high++;
            wins += low + (high - low + 1) / 2.0;
        }
        return wins / (targets.Length * (double)decoys.Length);
    }

    /// <summary>Half targets. The score sees one signal, and only the features carry the other.</summary>
    private static (List<double[]> Features, List<double> Scores, List<bool> IsDecoy, List<string> Keys) TwoSignals(int rows, int seed = 1)
    {
        var random = new Random(seed);
        var data = (Features: new List<double[]>(), Scores: new List<double>(), IsDecoy: new List<bool>(), Keys: new List<string>());
        for (int i = 0; i < rows; i++)
        {
            bool decoy = i % 2 == 1;
            double seen = (decoy ? 0 : 1) + Gaussian(random);
            double unseen = (decoy ? 0 : 1) + Gaussian(random);
            data.Features.Add([seen, unseen]);
            data.Scores.Add(seen);
            data.IsDecoy.Add(decoy);
            data.Keys.Add($"P{i}");
        }
        return data;
    }

    [Test]
    public void TheBlendAddsWhatTheScoreMisses()
    {
        var (features, scores, isDecoy, keys) = TwoSignals(3000);

        double[] blended = DiaTreeBlend.Blend(features, scores, isDecoy, keys, maxTrainingRows: 3000, seed: 0);

        Assert.That(Auc(blended, isDecoy), Is.GreaterThan(Auc(scores, isDecoy) + 0.03));
    }

    /// <summary>
    /// The blend decides the order and the score values stay the network's, so later stages see the same scale, and a
    /// decoy's score means what it did.
    /// </summary>
    [Test]
    public void TheBlendReordersTheScoresWithoutChangingTheirValues()
    {
        var (features, scores, isDecoy, keys) = TwoSignals(1000);

        double[] blended = DiaTreeBlend.Blend(features, scores, isDecoy, keys, maxTrainingRows: 1000, seed: 0);

        Assert.That(blended.Order(), Is.EqualTo(scores.Order()));
        Assert.That(blended, Is.Not.EqualTo(scores), "the order must actually change");
    }

    /// <summary>
    /// Whether a group is target or decoy is drawn at random, and the only feature names the group. A tree that saw a row's
    /// group while training could recite its label. Split by group, no tree sees it, so the ranking is chance.
    /// </summary>
    [Test]
    public void ARowIsNeverScoredByATreeThatSawItsGroup()
    {
        var random = new Random(3);
        var features = new List<double[]>();
        var scores = new List<double>();
        var isDecoy = new List<bool>();
        var keys = new List<string>();
        for (int group = 0; group < 1000; group++)
        {
            bool decoy = random.Next(2) == 0;
            for (int row = 0; row < 3; row++)
            {
                features.Add([group + 0.1 * random.NextDouble()]);
                scores.Add(random.NextDouble());
                isDecoy.Add(decoy);
                keys.Add($"G{group}");
            }
        }

        double[] blended = DiaTreeBlend.Blend(features, scores, isDecoy, keys, maxTrainingRows: 3000, seed: 0);

        Assert.That(Auc(blended, isDecoy), Is.EqualTo(0.5).Within(0.06));
    }

    /// <summary>
    /// Among the high scorers the feature points to targets, and among the many low scorers it points the other way. Trees
    /// taught on the top rows only rank the high scorers well; taught on every row, they learn the opposite.
    /// </summary>
    [Test]
    public void TheTreesLearnFromTheTopRowsByScore()
    {
        var random = new Random(5);
        var features = new List<double[]>();
        var scores = new List<double>();
        var isDecoy = new List<bool>();
        var keys = new List<string>();
        for (int i = 0; i < 4000; i++)
        {
            bool high = i < 1000;
            bool decoy = i % 2 == 1;
            bool featurePointsToTarget = high != decoy;
            features.Add([(featurePointsToTarget ? 1 : 0) + 0.5 * Gaussian(random)]);
            scores.Add((high ? 2 : 0) + random.NextDouble());
            isDecoy.Add(decoy);
            keys.Add($"P{i}");
        }
        var highRows = Enumerable.Range(0, 1000).ToList();
        double AucAmongHigh(double[] blended) => Auc(highRows.Select(i => blended[i]).ToList(), highRows.Select(i => isDecoy[i]).ToList());

        double topRows = AucAmongHigh(DiaTreeBlend.Blend(features, scores, isDecoy, keys, maxTrainingRows: 600, seed: 0));
        double everyRow = AucAmongHigh(DiaTreeBlend.Blend(features, scores, isDecoy, keys, maxTrainingRows: 4000, seed: 0));

        Assert.That(topRows, Is.GreaterThan(0.65));
        Assert.That(everyRow, Is.LessThan(0.5));
    }

    /// <summary>The A/B's identity checks need a seed to reproduce a run, however many threads it had.</summary>
    [Test]
    public void TheBlendIsTheSameForTheSameSeed()
    {
        var (features, scores, isDecoy, keys) = TwoSignals(1000);

        double[] first = DiaTreeBlend.Blend(features, scores, isDecoy, keys, maxTrainingRows: 1000, seed: 7, threads: 4);
        double[] second = DiaTreeBlend.Blend(features, scores, isDecoy, keys, maxTrainingRows: 1000, seed: 7, threads: 4);

        Assert.That(second, Is.EqualTo(first));
        Assert.That(first, Is.Not.EqualTo(scores));
    }

    /// <summary>A search too small to teach trees keeps the network's scores; it must not fail or reorder on noise.</summary>
    [TestCase(40, false, TestName = "TooFewRowsLeaveTheScoresAlone")]
    [TestCase(1000, true, TestName = "NoDecoysLeaveTheScoresAlone")]
    public void TooLittleToTrainOnLeavesTheScoresAlone(int rows, bool targetsOnly)
    {
        var (features, scores, isDecoy, keys) = TwoSignals(rows);
        if (targetsOnly)
            isDecoy = isDecoy.Select(_ => false).ToList();

        double[] blended = DiaTreeBlend.Blend(features, scores, isDecoy, keys, maxTrainingRows: rows, seed: 0);

        Assert.That(blended, Is.EqualTo(scores));
    }
}
