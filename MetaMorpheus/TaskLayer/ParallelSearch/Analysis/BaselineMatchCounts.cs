#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer;
using EngineLayer.SpectrumMatch;

namespace TaskLayer.ParallelSearch.Analysis;

/// <summary>Target and decoy counts at the confidence cutoffs over every PSM and every peptide of a search.</summary>
public readonly record struct GlobalMatchCounts(int PsmTargets, int PsmDecoys, int PeptideTargets, int PeptideDecoys);

/// <summary>
/// The baseline's PSM and peptide counts, and what is needed to update them for one transient search without sorting
/// and collapsing every PSM again.
/// <para>
/// A transient search changes only the scans it updated. The PSM counts are the baseline's, minus the updated scans'
/// baseline matches, plus their new ones. Peptides are collapsed per full sequence (<see cref="FilteredPsmsExtensions.CollapseToPeptides"/>),
/// so only sequences that an updated scan held before or holds now can change. Collapsing the PSMs of just those
/// sequences, in their original array order, gives exactly the representatives, and their order, that collapsing all
/// of them would. Transient peptides are taken from these, so they are peptides of scans this search changed. The full
/// lists also counted a baseline peptide as transient when its protein's accession was in the transient database
/// (a base database holding entries of a transient proteome, as the mixed regression base does), even on scans the
/// search never touched; that counted 3 transient peptides against 2 transient PSMs in one regression case.
/// </para>
/// <para>
/// Doing this once per run replaced about three sorts of every baseline PSM and a string grouping of all of them per
/// database (seconds and about 130 bytes per baseline PSM each time on deep files).
/// </para>
/// </summary>
public sealed class BaselineMatchCounts
{
    private const string NullSequenceKey = "\u0000";

    private readonly SpectralMatch?[] _baseline;
    private readonly Dictionary<string, List<int>> _indicesBySequence = new(StringComparer.Ordinal);
    private readonly Dictionary<string, (bool IsDecoy, bool Passes)> _baselinePeptides = new(StringComparer.Ordinal);
    private readonly GlobalMatchCounts _totals;

    /// <param name="baseline">The baseline PSM per scan (null where none), finished: FDR and parsimony done.</param>
    public BaselineMatchCounts(SpectralMatch?[] baseline, ConfidenceCutoff psmCutoff, ConfidenceCutoff peptideCutoff)
    {
        _baseline = baseline;

        int psmTargets = 0, psmDecoys = 0;
        for (int i = 0; i < baseline.Length; i++)
        {
            var psm = baseline[i];
            if (psm is null)
                continue;
            if (!_indicesBySequence.TryGetValue(Key(psm), out var indices))
                _indicesBySequence[Key(psm)] = indices = new List<int>();
            indices.Add(i);
            if (psmCutoff.Passes(psm))
            {
                if (psm.IsDecoy) psmDecoys++;
                else psmTargets++;
            }
        }

        int peptideTargets = 0, peptideDecoys = 0;
        foreach (var peptide in Collapse(baseline.Where(p => p is not null)!))
        {
            bool passes = peptideCutoff.Passes(peptide);
            _baselinePeptides[Key(peptide)] = (peptide.IsDecoy, passes);
            if (passes)
            {
                if (peptide.IsDecoy) peptideDecoys++;
                else peptideTargets++;
            }
        }

        _totals = new GlobalMatchCounts(psmTargets, psmDecoys, peptideTargets, peptideDecoys);
    }

    public GlobalMatchCounts BaselineTotals => _totals;

    /// <summary>
    /// The collapsed peptides of every full sequence that an updated scan held in the baseline or holds now, in the order
    /// collapsing all of <paramref name="allPsms"/> would list them.
    /// </summary>
    public List<SpectralMatch> CollapseAffectedPeptides(SpectralMatch?[] allPsms, IReadOnlyCollection<int> updatedIndexes,
        out HashSet<string> affectedSequences)
    {
        affectedSequences = new HashSet<string>(StringComparer.Ordinal);
        foreach (int i in updatedIndexes)
        {
            if (_baseline[i] is { } before)
                affectedSequences.Add(Key(before));
            if (allPsms[i] is { } after)
                affectedSequences.Add(Key(after));
        }

        var candidateIndexes = new HashSet<int>(updatedIndexes);
        foreach (string sequence in affectedSequences)
        {
            if (_indicesBySequence.TryGetValue(sequence, out var indices))
                candidateIndexes.UnionWith(indices);
        }

        var sequences = affectedSequences;
        var subset = candidateIndexes
            .Order()
            .Select(i => allPsms[i])
            .Where(p => p is not null && sequences.Contains(Key(p)))
            .Select(p => p!);

        return Collapse(subset).ToList();
    }

    /// <summary>
    /// The run-wide counts for one transient search, read at the cutoffs now, after parsimony has trimmed the transient
    /// matches as the full lists would have been read.
    /// </summary>
    public GlobalMatchCounts Count(SpectralMatch?[] allPsms, IReadOnlyCollection<int> updatedIndexes,
        List<SpectralMatch> affectedPeptides, HashSet<string> affectedSequences,
        ConfidenceCutoff psmCutoff, ConfidenceCutoff peptideCutoff)
    {
        int psmTargets = _totals.PsmTargets, psmDecoys = _totals.PsmDecoys;
        foreach (int i in updatedIndexes)
        {
            Tally(_baseline[i], psmCutoff, -1, ref psmTargets, ref psmDecoys);
            Tally(allPsms[i], psmCutoff, +1, ref psmTargets, ref psmDecoys);
        }

        int peptideTargets = _totals.PeptideTargets, peptideDecoys = _totals.PeptideDecoys;
        foreach (string sequence in affectedSequences)
        {
            if (_baselinePeptides.TryGetValue(sequence, out var before) && before.Passes)
            {
                if (before.IsDecoy) peptideDecoys--;
                else peptideTargets--;
            }
        }
        foreach (var peptide in affectedPeptides)
            Tally(peptide, peptideCutoff, +1, ref peptideTargets, ref peptideDecoys);

        return new GlobalMatchCounts(psmTargets, psmDecoys, peptideTargets, peptideDecoys);
    }

    private static void Tally(SpectralMatch? match, ConfidenceCutoff cutoff, int sign, ref int targets, ref int decoys)
    {
        if (match is null || !cutoff.Passes(match))
            return;
        if (match.IsDecoy) decoys += sign;
        else targets += sign;
    }

    /// <summary>As the full path: sort by the default comparer (stable, so array order breaks ties), then collapse.</summary>
    private static IEnumerable<SpectralMatch> Collapse(IEnumerable<SpectralMatch> psmsInArrayOrder) =>
        psmsInArrayOrder.OrderByDescending(p => p).CollapseToPeptides(true);

    private static string Key(SpectralMatch match) => match.FullSequence ?? NullSequenceKey;
}
