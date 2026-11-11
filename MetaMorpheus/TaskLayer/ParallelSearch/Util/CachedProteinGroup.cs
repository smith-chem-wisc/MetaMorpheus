using EngineLayer;
using Omics;
using System.Collections.Generic;

namespace TaskLayer.ParallelSearch;

/// <summary>
/// Represents a cached version of a ProteinGroup for use in parallel search.
/// </summary>
internal class CachedProteinGroup(ProteinGroup proteinGroup)
{
    private ProteinGroup ProteinGroup { get; } = BuildCached(proteinGroup);

    private static ProteinGroup BuildCached(ProteinGroup pg)
    {
        var cached = new ProteinGroup(
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
        cached.SequenceCoverageFraction.AddRange(pg.SequenceCoverageFraction);
        cached.SequenceCoverageDisplayList.AddRange(pg.SequenceCoverageDisplayList);
        cached.SequenceCoverageDisplayListWithMods.AddRange(pg.SequenceCoverageDisplayListWithMods);
        cached.FragmentSequenceCoverageDisplayList.AddRange(pg.FragmentSequenceCoverageDisplayList);
        return cached;
    }

    public ProteinGroup GetProteinGroup() => ProteinGroup;

    public ProteinGroup CreateRuntimeCopy()
    {
        var copy = new ProteinGroup(
            [..ProteinGroup.Proteins],
            [..ProteinGroup.AllPeptides],
            [..ProteinGroup.UniquePeptides])
        {
            AllPsmsBelowOnePercentFDR = [..ProteinGroup.AllPsmsBelowOnePercentFDR],
            ProteinGroupScore = ProteinGroup.ProteinGroupScore,
            BestPeptideScore = ProteinGroup.BestPeptideScore,
            BestPeptideQValue = ProteinGroup.BestPeptideQValue,
            BestPeptidePEP = ProteinGroup.BestPeptidePEP,
            QValue = ProteinGroup.QValue,
            CumulativeTarget = ProteinGroup.CumulativeTarget,
            CumulativeDecoy = ProteinGroup.CumulativeDecoy,
            DisplayModsOnPeptides = ProteinGroup.DisplayModsOnPeptides,
            // SampleGroupResults (quant + occupancy) is derived from this and rebuilt on read, so it is not copied
            OccupancySitesToReport = ProteinGroup.OccupancySitesToReport,
        };
        copy.SequenceCoverageFraction.AddRange(ProteinGroup.SequenceCoverageFraction);
        copy.SequenceCoverageDisplayList.AddRange(ProteinGroup.SequenceCoverageDisplayList);
        copy.SequenceCoverageDisplayListWithMods.AddRange(ProteinGroup.SequenceCoverageDisplayListWithMods);
        copy.FragmentSequenceCoverageDisplayList.AddRange(ProteinGroup.FragmentSequenceCoverageDisplayList);
        return copy;
    }
}
