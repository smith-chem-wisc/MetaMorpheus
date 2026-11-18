#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// Protein-level support for a precursor: how many other peptides of its protein are confident in the run. A peptide of a
/// protein whose other peptides are found is more likely there, and a lone peptide of a protein nothing else supports less.
/// DIA-NN 2.3.2 re-scores on protein information after its network.
/// </summary>
/// <remarks>
/// A decoy counts toward its target's protein (its accession less <see cref="DecoyPrefix"/>), so the decoys that set q get the
/// support a false hit in that protein would. An entrapment benchmark must do the same for its entrapments, or support
/// removes them for free: an entrapment in a protein of its own is never in a present protein, while a real false hit falls
/// in a random one (see <see cref="DiaLibrarySearchParameters.ProteinKey"/>).
/// </remarks>
public static class DiaProteinSupport
{
    /// <summary>The prefix of a decoy protein's accession.</summary>
    public const string DecoyPrefix = "DECOY_";

    /// <summary>One precursor as support reads it.</summary>
    /// <param name="Accession">Its protein accessions, '|' or ';' joined; empty or null for none.</param>
    /// <param name="BaseSequence">Its peptide: charge states and modified forms of one peptide count once.</param>
    /// <param name="IsDecoy">A decoy is never confident and counts toward its target's protein.</param>
    /// <param name="Confident">A target-side precursor that passes the provisional threshold.</param>
    public sealed record Row(string? Accession, string BaseSequence, bool IsDecoy, bool Confident);

    /// <summary>
    /// For each row, the most, over its proteins, of the confident target-side peptides of that protein other than its own.
    /// </summary>
    /// <param name="key">The protein an accession counts toward; null takes the accession as it is.</param>
    /// <exception cref="ArgumentNullException"><paramref name="rows"/> is null.</exception>
    public static double[] Measure(IReadOnlyList<Row> rows, Func<string, string>? key = null)
    {
        ArgumentNullException.ThrowIfNull(rows);
        var proteins = new string[rows.Count][];
        var confidentPeptides = new Dictionary<string, HashSet<string>>(StringComparer.Ordinal);
        for (int i = 0; i < rows.Count; i++)
        {
            proteins[i] = ProteinsOf(rows[i].Accession, key);
            if (rows[i].IsDecoy || !rows[i].Confident)
                continue;
            foreach (string protein in proteins[i])
            {
                if (!confidentPeptides.TryGetValue(protein, out var peptides))
                    confidentPeptides[protein] = peptides = new HashSet<string>(StringComparer.Ordinal);
                peptides.Add(rows[i].BaseSequence);
            }
        }
        var support = new double[rows.Count];
        for (int i = 0; i < rows.Count; i++)
        {
            foreach (string protein in proteins[i])
            {
                if (!confidentPeptides.TryGetValue(protein, out var peptides))
                    continue;
                int others = peptides.Count - (!rows[i].IsDecoy && peptides.Contains(rows[i].BaseSequence) ? 1 : 0);
                support[i] = Math.Max(support[i], others);
            }
        }
        return support;
    }

    /// <summary>The proteins an accession string names: split on '|' and ';', a decoy's prefix dropped, each through <paramref name="key"/>.</summary>
    private static string[] ProteinsOf(string? accession, Func<string, string>? key)
    {
        if (string.IsNullOrEmpty(accession))
            return [];
        return accession.Split(['|', ';'], StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
            .Select(a => a.StartsWith(DecoyPrefix, StringComparison.Ordinal) ? a[DecoyPrefix.Length..] : a)
            .Where(a => a.Length > 0)
            .Select(a => key is null ? a : key(a))
            .Distinct(StringComparer.Ordinal).ToArray();
    }
}
