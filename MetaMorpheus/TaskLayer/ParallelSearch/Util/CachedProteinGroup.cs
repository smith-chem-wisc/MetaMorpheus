using EngineLayer;
using Omics;
using System.Collections.Generic;

namespace TaskLayer.ParallelSearch;

/// <summary>
/// Represents a cached version of a ProteinGroup for use in parallel search.
/// </summary>
internal class CachedProteinGroup(ProteinGroup proteinGroup)
{
    private ProteinGroup ProteinGroup { get; } = DeepCopy(proteinGroup);

    /// <summary>
    /// Copies a group so that writing it, which changes the group, leaves the original untouched.
    /// </summary>
    public static ProteinGroup DeepCopy(ProteinGroup pg)
    {
        var copy = new ProteinGroup(
            [..pg.Proteins],
            [..pg.AllPeptides],
            [..pg.UniquePeptides])
        {
            AllPsmsBelowOnePercentFDR = [..pg.AllPsmsBelowOnePercentFDR],
            ProteinGroupScore = pg.ProteinGroupScore,
            BestPeptideScore = pg.BestPeptideScore,
            BestPeptideQValue = pg.BestPeptideQValue,
            BestPeptidePEP = pg.BestPeptidePEP,
            QValue = pg.QValue,
            CumulativeTarget = pg.CumulativeTarget,
            CumulativeDecoy = pg.CumulativeDecoy,
            DisplayModsOnPeptides = pg.DisplayModsOnPeptides,
            // SampleGroupResults (quant + occupancy) is derived from this and rebuilt on read, so it is not copied
            OccupancySitesToReport = pg.OccupancySitesToReport,
        };
        copy.SequenceCoverageFraction.AddRange(pg.SequenceCoverageFraction);
        copy.SequenceCoverageDisplayList.AddRange(pg.SequenceCoverageDisplayList);
        copy.SequenceCoverageDisplayListWithMods.AddRange(pg.SequenceCoverageDisplayListWithMods);
        copy.FragmentSequenceCoverageDisplayList.AddRange(pg.FragmentSequenceCoverageDisplayList);
        return copy;
    }

    public ProteinGroup GetProteinGroup() => ProteinGroup;

    public ProteinGroup CreateRuntimeCopy() => DeepCopy(ProteinGroup);
}
