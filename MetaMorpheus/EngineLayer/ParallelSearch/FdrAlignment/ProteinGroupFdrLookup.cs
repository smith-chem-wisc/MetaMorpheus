#nullable enable

using System.Collections.Generic;
using EngineLayer.SpectrumMatch;

namespace EngineLayer.ParallelSearch.FdrAlignment;

public sealed class ProteinGroupFdrLookup : BaselineRankLookup<ProteinGroupRank, ProteinGroupFdr>
{
    public ProteinGroupFdrLookup(IEnumerable<(ProteinGroupRank Rank, ProteinGroupFdr Fdr)> baseline, FilterType filterType)
        : base(baseline, ProteinGroupRank.BetterFirst)
    {
        FilterType = filterType;
    }

    /// <summary>
    /// The filter the baseline's protein FDR ranked groups by. Transient groups must be ranked the same way.
    /// </summary>
    public FilterType FilterType { get; }
}
