#nullable enable
using System.Collections.Generic;
using Chromatography.RetentionTimeCalibration;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using Omics;

namespace TaskLayer.ParallelSearch.Analysis;

/// <summary>
/// Context containing all data needed for analysis of a single transient database
/// Provides read-only access to search results and metadata
/// </summary>
public class TransientDatabaseContext : ITransientDatabaseResult
{
    public string DatabaseName { get; init; } = string.Empty;
    public DbForTask TransientDatabase { get; init; } = null!;
    public List<IBioPolymer> TransientProteins { get; init; } = null!;
    public HashSet<string> TransientProteinAccessions { get; init; } = null!;

    // PSM and Peptide data
    public List<SpectralMatch> AllPsms { get; init; } = null!;
    public List<SpectralMatch> TransientPsms { get; init; } = null!;
    public List<SpectralMatch> AllPeptides { get; init; } = null!;
    public List<SpectralMatch> TransientPeptides { get; init; } = null!;

    // Protein data (may be null if parsimony not performed)
    public List<ProteinGroup>? ProteinGroups { get; init; }
    public List<ProteinGroup>? TransientProteinGroups { get; init; }

    // Chronologer % acetonitrile to minutes, per spectra file path, fitted once on the baseline search; null when not fitted
    public IReadOnlyDictionary<string, IrtCalibrationModel>? RetentionTimeCalibrations { get; init; }

    // Parameters and metadata
    public CommonParameters CommonParameters { get; init; } = null!;
    public int TotalProteins { get; init; }
    public int TransientPeptideCount { get; init; }

    public string OutputFolder { get; init; } = string.Empty;
    public List<string> NestedIds { get; init; } = [];

    // the confidence cutoffs a normal search would use for these results, decided once and shared by every collector
    private ConfidenceCutoff? _psmCutoff;
    private ConfidenceCutoff? _peptideCutoff;
    private double? _proteinGroupQValueThreshold;
    public ConfidenceCutoff PsmCutoff => _psmCutoff ??= ConfidenceCutoff.For(AllPsms, CommonParameters, peptideLevel: false);
    public ConfidenceCutoff PeptideCutoff => _peptideCutoff ??= ConfidenceCutoff.For(AllPeptides, CommonParameters, peptideLevel: true);
    public double ProteinGroupQValueThreshold =>
        _proteinGroupQValueThreshold ??= ConfidenceCutoff.ProteinGroupThreshold(AllPsms, CommonParameters);
}
