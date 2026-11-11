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
    public const string ProteinGroupBacterialTargets = "ProteinGroupBacterialTargets";
    public const string ProteinGroupBacterialDecoys = "ProteinGroupBacterialDecoys";
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
        yield return ProteinGroupBacterialTargets;
        yield return ProteinGroupBacterialDecoys;
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
        double qValueThreshold = Math.Min(context.CommonParameters.QValueThreshold, context.CommonParameters.PepQValueThreshold);
        var totalTargets = context.ProteinGroups!.Count(p => !p.IsDecoy);
        var totalDecoys = context.ProteinGroups!.Count(p => p.IsDecoy);
        var totalTransientTargets = context.TransientProteinGroups!.Count(p => !p.IsDecoy);


        var (targetsGlobalAtCutoff, decoysGlobalAtCutoff, unambiguousTargetsGlobal, unambiguousDecoysGlobal) =
            AnalyzeProteinGroups(context.ProteinGroups!, qValueThreshold);

        var (bacterialTargetsAtCutoff, bacterialDecoysAtCutoff, unambiguousTargets, unambiguousDecoys) =
            AnalyzeProteinGroups(context.TransientProteinGroups!, qValueThreshold);

        var transientTargets = context.TransientProteinGroups!.Where(p => !p.IsDecoy).ToList();

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
            [ProteinGroupBacterialTargets] = totalTransientTargets,
            [ProteinGroupBacterialDecoys] = context.TransientProteinGroups!.Count - totalTransientTargets,
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

    private (int Targets, int Decoys, int UnambiguousTargets, int UnambiguousDecoys)
        AnalyzeProteinGroups(List<ProteinGroup> proteinGroups, double qValueThreshold)
    {
        int targets = 0, decoys = 0, unambiguousTargets = 0, unambiguousDecoys = 0;

        foreach (var pg in proteinGroups)
        {
            if (pg.QValue > qValueThreshold || pg.AllPeptides.Count < 2)
                continue;

            bool isOrganismAmbiguous = false;
            foreach (var protein in pg.Proteins)
            {
                string? organism = protein.Organism;
                if (organism != null && organism.Contains(_targetOrganism, StringComparison.Ordinal))
                {
                    isOrganismAmbiguous = true;
                    break;
                }
            }

            if (pg.IsDecoy)
            {
                decoys++;
                if (!isOrganismAmbiguous)
                    unambiguousDecoys++;
            }
            else
            {
                targets++;
                if (!isOrganismAmbiguous)
                    unambiguousTargets++;
            }
        }

        return (targets, decoys, unambiguousTargets, unambiguousDecoys);
    }
}
