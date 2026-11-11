#nullable enable
using System.Collections.Generic;
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

    // Parameters and metadata
    public CommonParameters CommonParameters { get; init; } = null!;
    public int TotalProteins { get; init; }
    public int TransientPeptideCount { get; init; }

    public string OutputFolder { get; init; } = string.Empty;
    public List<string> NestedIds { get; init; } = [];
}
