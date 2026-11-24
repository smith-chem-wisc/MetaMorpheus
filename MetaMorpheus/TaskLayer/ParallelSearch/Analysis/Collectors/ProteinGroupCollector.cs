#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using MathNet.Numerics.Statistics;

namespace TaskLayer.ParallelSearch.Analysis.Collectors;

/// <summary>
/// Analyzer for protein group metrics
/// Only runs if parsimony was performed
/// </summary>
public class ProteinGroupCollector : IMetricCollector
{
    // Output Column Names
    public const string TargetProteinGroupsAtQValueThreshold = "TargetProteinGroupsAtQValueThreshold";
    public const string TargetProteinGroupsFromTransientDb = "TargetProteinGroupsFromTransientDb";
    public const string TargetProteinGroupsFromTransientDbAtQValueThreshold = "TargetProteinGroupsFromTransientDbAtQValueThreshold";
    public const string DecoyProteinGroupsFromTransientDbAtQValueThreshold = "DecoyProteinGroupsFromTransientDbAtQValueThreshold";
    public const string ProteinGroupBacterialUnambiguousTargets = "ProteinGroupBacterialUnambiguousTargets";
    public const string ProteinGroupBacterialUnambiguousDecoys = "ProteinGroupBacterialUnambiguousDecoys";
    public const string ProteinGroupTargets = "ProteinGroupTargets";
    public const string ProteinGroupDecoys = "ProteinGroupDecoys";
    public const string MedianPsmsPerProteinGroup = "MedianPsmsPerProteinGroup";
    public const string MedianPeptidesPerProteinGroup = "MedianPeptidesPerProteinGroup";
    public const string MedianUniquePeptidesPerProteinGroup = "MedianUniquePeptidesPerProteinGroup";
    public const string AllPeptidesPerProteinGroup = "AllPeptidesPerProteinGroup";
    public const string AllUniquePeptidesPerProteinGroup = "AllUniquePeptidesPerProteinGroup";
    public const string AllPsmsPerProteinGroup = "AllPsmsPerProteinGroup";
    public const string AllSequenceCoverageFractions = "AllSequenceCoverageFractions";
    public const string AllFragmentSequenceCoverageFractions = "AllFragmentSequenceCoverageFractions";

    private readonly string _targetOrganism;

    public ProteinGroupCollector(string targetOrganism = "Homo sapiens")
    {
        _targetOrganism = targetOrganism;
    }

    public string CollectorName => "ProteinGroups";

    public IEnumerable<string> GetOutputColumns()
    {
        yield return ProteinGroupTargets;
        yield return ProteinGroupDecoys;
        yield return TargetProteinGroupsAtQValueThreshold;
        yield return TargetProteinGroupsFromTransientDb;
        yield return TargetProteinGroupsFromTransientDbAtQValueThreshold;
        yield return DecoyProteinGroupsFromTransientDbAtQValueThreshold;
        yield return ProteinGroupBacterialUnambiguousTargets;
        yield return ProteinGroupBacterialUnambiguousDecoys;
        yield return MedianPeptidesPerProteinGroup;
        yield return AllPeptidesPerProteinGroup;
        yield return MedianUniquePeptidesPerProteinGroup;
        yield return AllUniquePeptidesPerProteinGroup;
        yield return MedianPsmsPerProteinGroup;
        yield return AllPsmsPerProteinGroup;
        yield return AllSequenceCoverageFractions;
        yield return AllFragmentSequenceCoverageFractions;
    }

    public bool CanCollectData(TransientDatabaseContext context)
    {
        return context.ProteinGroups != null && context.TransientProteinGroups != null;
    }

    public Dictionary<string, object> CollectData(TransientDatabaseContext context)
    {
        double qValueThreshold = context.ProteinGroupQValueThreshold;
        var totalTargets = context.ProteinGroups!.Count(p => !p.IsDecoy);
        var totalDecoys = context.ProteinGroups!.Count(p => p.IsDecoy);
        var totalTransientTargets = context.TransientProteinGroups!.Count(p => !p.IsDecoy);


        var (targetsGlobalAtCutoff, decoysGlobalAtCutoff, unambiguousTargetsGlobal, unambiguousDecoysGlobal) =
            AnalyzeProteinGroups(context.ProteinGroups!, qValueThreshold, IsForegroundOrganism);

        // a transient group is ambiguous when it holds a foreground-organism protein or any protein outside the database
        Func<ProteinGroup, bool> isAmbiguous = pg => IsForegroundOrganism(pg) || context.IsSharedWithBaseline(pg);
        var (bacterialTargetsAtCutoff, bacterialDecoysAtCutoff, unambiguousTargets, unambiguousDecoys) =
            AnalyzeProteinGroups(context.TransientProteinGroups!, qValueThreshold, isAmbiguous);

        // the per-group distributions describe the groups counted as unambiguous targets above: confident groups that are
        // this database's own (unconfident groups gave databases with no confident PSMs KS p-values of 1e-5)
        var transientTargets = context.TransientProteinGroups!
            .Where(pg => !pg.IsDecoy && IsConfident(pg, qValueThreshold) && !isAmbiguous(pg))
            .ToList();

        var peptidesPerProtein = transientTargets.Select(pg => (double)pg.AllPeptides.Count(p => !p.Parent.IsDecoy));
        var uniquePeptidesPerProtein = transientTargets.Select(pg => (double)pg.UniquePeptides.Count(p => !p.Parent.IsDecoy));
        var psmsPerProtein = transientTargets.Select(pg => (double)pg.AllPsmsBelowOnePercentFDR.Count(p => !p.IsDecoy));

        var seqCoverage = transientTargets.SelectMany(pg => pg.SequenceCoverageFraction.Where(v => v > 0));
        var fragSeqCoverage = transientTargets.SelectMany(pg => pg.FragmentSequenceCoverageDisplayList
            .Select(s => s.Length > 0 ? s.Count(char.IsUpper) / (double)s.Length : 0.0));

        return new Dictionary<string, object>
        {
            [ProteinGroupTargets] = totalTargets,
            [ProteinGroupDecoys] = totalDecoys,
            [TargetProteinGroupsAtQValueThreshold] = targetsGlobalAtCutoff,
            [TargetProteinGroupsFromTransientDb] = totalTransientTargets,
            [TargetProteinGroupsFromTransientDbAtQValueThreshold] = bacterialTargetsAtCutoff,
            [DecoyProteinGroupsFromTransientDbAtQValueThreshold] = bacterialDecoysAtCutoff,
            [ProteinGroupBacterialUnambiguousTargets] = unambiguousTargets,
            [ProteinGroupBacterialUnambiguousDecoys] = unambiguousDecoys,

            [MedianPeptidesPerProteinGroup] = peptidesPerProtein.Median(),
            [AllPeptidesPerProteinGroup] = peptidesPerProtein.ToArray(),
            [MedianUniquePeptidesPerProteinGroup] = uniquePeptidesPerProtein.Median(),
            [AllUniquePeptidesPerProteinGroup] = uniquePeptidesPerProtein.ToArray(),
            [MedianPsmsPerProteinGroup] = psmsPerProtein.Median(),
            [AllPsmsPerProteinGroup] = psmsPerProtein.ToArray(),

            [AllSequenceCoverageFractions] = seqCoverage.ToArray(),
            [AllFragmentSequenceCoverageFractions] = fragSeqCoverage.ToArray(),
        };
    }

    private static bool IsConfident(ProteinGroup pg, double qValueThreshold) =>
        pg.QValue <= qValueThreshold && pg.AllPeptides.Count >= 2;

    private bool IsForegroundOrganism(ProteinGroup pg) =>
        pg.Proteins.Any(p => p.Organism is { } organism && organism.Contains(_targetOrganism, StringComparison.Ordinal));

    private (int Targets, int Decoys, int UnambiguousTargets, int UnambiguousDecoys)
        AnalyzeProteinGroups(List<ProteinGroup> proteinGroups, double qValueThreshold, Func<ProteinGroup, bool> isAmbiguous)
    {
        int targets = 0, decoys = 0, unambiguousTargets = 0, unambiguousDecoys = 0;

        foreach (var pg in proteinGroups)
        {
            if (!IsConfident(pg, qValueThreshold))
                continue;

            bool isAmbiguousGroup = isAmbiguous(pg);

            if (pg.IsDecoy)
            {
                decoys++;
                if (!isAmbiguousGroup)
                    unambiguousDecoys++;
            }
            else
            {
                targets++;
                if (!isAmbiguousGroup)
                    unambiguousTargets++;
            }
        }

        return (targets, decoys, unambiguousTargets, unambiguousDecoys);
    }
}
