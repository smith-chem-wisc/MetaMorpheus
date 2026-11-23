#nullable enable
using System.Collections.Generic;
using System.Linq;
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

    /// <summary>The matches that tied with a baseline protein before transient parsimony trimmed their hypotheses.</summary>
    public BaselineTies BaselineTies { get; init; } = BaselineTies.Empty;

    /// <summary>
    /// Run-wide target and decoy counts at the cutoffs, when they were updated from the baseline's instead of counted
    /// over <see cref="AllPsms"/> and <see cref="AllPeptides"/>; those lists are then empty unless an output needs them.
    /// Read them through <see cref="GetGlobalCounts"/>.
    /// </summary>
    public GlobalMatchCounts? GlobalCounts { get; init; }

    /// <summary>The run-wide counts: <see cref="GlobalCounts"/> when set, otherwise counted over the full lists.</summary>
    public GlobalMatchCounts GetGlobalCounts() => GlobalCounts ?? new GlobalMatchCounts(
        AllPsms.Count(p => !p.IsDecoy && PsmCutoff.Passes(p)),
        AllPsms.Count(p => p.IsDecoy && PsmCutoff.Passes(p)),
        AllPeptides.Count(p => !p.IsDecoy && PeptideCutoff.Passes(p)),
        AllPeptides.Count(p => p.IsDecoy && PeptideCutoff.Passes(p)));

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

    /// <summary>
    /// True when any of the match's hypotheses comes from a protein outside this transient database (the baseline, its
    /// contaminants, or another organism in it): the spectrum is explained without this database, so it is no evidence
    /// for it. I and L have the same mass, so a transient peptide that is an I/L isomer of a baseline peptide ties with it
    /// and is caught here too (on PXD038928, one phage peptide matched 748 human PSMs this way). Ties that parsimony
    /// removed from the match are read from <see cref="BaselineTies"/>.
    /// </summary>
    public bool IsSharedWithBaseline(SpectralMatch match) =>
        TransientProteinAccessions is not null
        && (match.HypothesesInAddedOrder.Any(h => !TransientProteinAccessions.Contains(h.SpecificBioPolymer.Parent.Accession))
            || BaselineTies.Contains(match, TransientProteinAccessions));

    /// <summary>
    /// True when any protein in the group comes from outside this transient database, or when none of the group's PSMs is
    /// the database's own. Parsimony can put a transient protein in a group of its own while every one of its PSMs ties
    /// with a baseline protein (HSV-1 UP000180758 against the strain-17 entries in the mixed regression base: 0 unambiguous
    /// peptides, 5 groups).
    /// </summary>
    public bool IsSharedWithBaseline(ProteinGroup group) =>
        TransientProteinAccessions is not null
        && (group.Proteins.Any(p => !TransientProteinAccessions.Contains(p.Accession))
            || !group.AllPsmsBelowOnePercentFDR.OfType<SpectralMatch>().Any(psm => !IsSharedWithBaseline(psm)));

    private List<SpectralMatch>? _confidentUnambiguousPsms;
    private List<SpectralMatch>? _confidentUnambiguousPeptides;

    /// <summary>Transient PSMs (targets and decoys) that pass the PSM cutoff and are not shared with the baseline.</summary>
    public List<SpectralMatch> ConfidentUnambiguousPsms => _confidentUnambiguousPsms ??=
        TransientPsms.Where(p => PsmCutoff.Passes(p) && !IsSharedWithBaseline(p)).ToList();

    /// <summary>Transient peptides (targets and decoys) that pass the peptide cutoff and are not shared with the baseline.</summary>
    public List<SpectralMatch> ConfidentUnambiguousPeptides => _confidentUnambiguousPeptides ??=
        TransientPeptides.Where(p => PeptideCutoff.Passes(p) && !IsSharedWithBaseline(p)).ToList();
}
