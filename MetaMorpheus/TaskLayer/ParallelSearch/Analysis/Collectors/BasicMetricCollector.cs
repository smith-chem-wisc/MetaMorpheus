#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

namespace TaskLayer.ParallelSearch.Analysis.Collectors;

/// <summary>
/// Analyzer for basic search metrics (PSM counts, protein counts, etc.)
/// This replaces the hardcoded metrics currently in TransientDatabaseSearchResults
/// </summary>
public class BasicMetricCollector : IMetricCollector
{
    // Column Names
    public const string TotalProteins = "TotalProteins";
    public const string TransientProteinCount = "TransientProteinCount";
    public const string TransientPeptideCount = "TransientPeptideCount";
    public const string TargetPsmsAtQValueThreshold = "TargetPsmsAtQValueThreshold";
    public const string TargetPsmsFromTransientDb = "TargetPsmsFromTransientDb";
    public const string TargetPsmsFromTransientDbAtQValueThreshold = "TargetPsmsFromTransientDbAtQValueThreshold";
    public const string TargetPeptidesAtQValueThreshold = "TargetPeptidesAtQValueThreshold";
    public const string TargetPeptidesFromTransientDb = "TargetPeptidesFromTransientDb";
    public const string TargetPeptidesFromTransientDbAtQValueThreshold = "TargetPeptidesFromTransientDbAtQValueThreshold";

    public string CollectorName => "ResultCount";

    public IEnumerable<string> GetOutputColumns()
    {
        yield return TotalProteins;
        yield return TransientProteinCount;
        yield return TransientPeptideCount;
        yield return TargetPsmsAtQValueThreshold;
        yield return TargetPsmsFromTransientDb;
        yield return TargetPsmsFromTransientDbAtQValueThreshold;
        yield return TargetPeptidesAtQValueThreshold;
        yield return TargetPeptidesFromTransientDb;
        yield return TargetPeptidesFromTransientDbAtQValueThreshold;
    }

    public bool CanCollectData(TransientDatabaseContext context)
    {
        return context.AllPsms != null
            && context.TransientPsms != null
            && context.AllPeptides != null
            && context.TransientPeptides != null;
    }

    public Dictionary<string, object> CollectData(TransientDatabaseContext context)
    {
        double qValueThreshold = Math.Min(context.CommonParameters.QValueThreshold, context.CommonParameters.PepQValueThreshold);

        return new Dictionary<string, object>
        {
            [TotalProteins] = context.TotalProteins,
            [TransientProteinCount] = context.TransientProteinAccessions.Count,
            [TransientPeptideCount] = context.TransientPeptideCount,

            [TargetPsmsAtQValueThreshold] = context.AllPsms.Count(p => !p.IsDecoy && p.GetFdrInfo(false)?.QValue <= qValueThreshold),
            [TargetPsmsFromTransientDb] = context.TransientPsms.Count(p => !p.IsDecoy),
            [TargetPsmsFromTransientDbAtQValueThreshold] = context.TransientPsms.Count(p => !p.IsDecoy && p.GetFdrInfo(false)?.QValue <= qValueThreshold),

            [TargetPeptidesAtQValueThreshold] = context.AllPeptides.Count(p => !p.IsDecoy && p.PeptideFdrInfo?.QValue <= qValueThreshold),
            [TargetPeptidesFromTransientDb] = context.TransientPeptides.Count(p => !p.IsDecoy),
            [TargetPeptidesFromTransientDbAtQValueThreshold] = context.TransientPeptides.Count(p => !p.IsDecoy && p.PeptideFdrInfo?.QValue <= qValueThreshold)
        };
    }
}