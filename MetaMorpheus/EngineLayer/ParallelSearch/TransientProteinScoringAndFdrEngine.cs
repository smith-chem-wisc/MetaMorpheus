using System;
using System.Collections.Generic;
using System.Linq;
using Easy.Common.Extensions;
using EngineLayer.ParallelSearch.FdrAlignment;
using MzLibUtil;

namespace EngineLayer.ParallelSearch;

/// <summary>
/// Protein Scoring engine that only scores the transient protein groups, but applies FDR to the transient and baseline protein groups together. 
/// This is used for parallel searches where the baseline search has already been scored and had FDR applied, and we want to score the transient search and apply FDR to the combined set of baseline and transient protein groups.
/// </summary>
public sealed class TransientProteinScoringAndFdrEngine : ProteinScoringAndFdrEngine
{
    private readonly List<ProteinGroup> _baselineProteinGroups;
    private readonly List<ProteinGroup> _transientProteinGroups;
    private readonly IEnumerable<SpectralMatch> _scoringPsms;
    private readonly ProteinGroupFdrAlignmentService _proteinGroupFdrAlignmentService;

    public TransientProteinScoringAndFdrEngine(
        List<ProteinGroup> baselineProteinGroups,
        List<ProteinGroup> transientProteinGroups,
        List<SpectralMatch> neighborhoodPsms,
        ProteinGroupFdrAlignmentService proteinGroupFdrAlignmentService,
        bool noOneHitWonders,
        bool treatModPeptidesAsDifferentPeptides,
        bool mergeIndistinguishableProteinGroups,
        CommonParameters commonParameters,
        List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
        List<string> nestedIds)
        : base(
            transientProteinGroups,
            neighborhoodPsms,
            psmsAlreadyFiltered: true,
            noOneHitWonders,
            treatModPeptidesAsDifferentPeptides,
            mergeIndistinguishableProteinGroups,
            commonParameters,
            fileSpecificParameters,
            nestedIds)
    {
        ArgumentNullException.ThrowIfNull(baselineProteinGroups);
        ArgumentNullException.ThrowIfNull(transientProteinGroups);
        ArgumentNullException.ThrowIfNull(neighborhoodPsms);
        ArgumentNullException.ThrowIfNull(proteinGroupFdrAlignmentService);

        _baselineProteinGroups = baselineProteinGroups;
        _transientProteinGroups = transientProteinGroups;
        _scoringPsms = neighborhoodPsms;
        _proteinGroupFdrAlignmentService = proteinGroupFdrAlignmentService;
    }

    protected override MetaMorpheusEngineResults RunSpecific()
    {
        ProteinScoringAndFdrResults analysisResults = new(this);

        if (_transientProteinGroups.Count > 0)
        {
            ScoreProteinGroups(_transientProteinGroups, _scoringPsms);
        }

        List<ProteinGroup> filteredTransientGroups = ApplyNoOneHitWondersFilter(_transientProteinGroups);
        PopulateBestPeptideMetrics(filteredTransientGroups);

        if (_proteinGroupFdrAlignmentService.HasBaselineCache)
        {
            _ = _proteinGroupFdrAlignmentService.ApplyBaseline(filteredTransientGroups);
        }
        else if (filteredTransientGroups.Count > 0)
        {
            filteredTransientGroups = DoProteinFdr(filteredTransientGroups);
        }

        analysisResults.SortedAndScoredProteinGroups = _baselineProteinGroups.Concat(filteredTransientGroups)
            .OrderBy(pg => pg.QValue)
            .ThenByDescending(pg => pg.ProteinGroupScore)
            .ToList();

        // Calculate sequence coverage only for those results who need it. 
        foreach (var proteinGroup in analysisResults.SortedAndScoredProteinGroups.Where(p => p.SequenceCoverageFraction.IsNullOrEmpty()))
            proteinGroup.CalculateSequenceCoverage();

        return analysisResults;
    }
}
