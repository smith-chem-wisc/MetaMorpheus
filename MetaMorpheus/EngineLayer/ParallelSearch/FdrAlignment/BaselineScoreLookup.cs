#nullable enable

using System.Collections.Generic;

namespace EngineLayer.ParallelSearch.FdrAlignment;

/// <summary>
/// PSM and peptide FDR rank by score, highest first.
/// </summary>
public sealed class BaselineScoreLookup<T> : BaselineRankLookup<double, T>
{
    private static readonly IComparer<double> HigherScoreFirst = Comparer<double>.Create((a, b) => b.CompareTo(a));

    public BaselineScoreLookup(IEnumerable<(double Score, T Value)> baseline)
        : base(baseline, HigherScoreFirst)
    {
    }
}
