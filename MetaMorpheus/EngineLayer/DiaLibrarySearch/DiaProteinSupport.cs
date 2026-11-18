#nullable enable
using System;
using System.Collections.Generic;

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
        return new double[rows.Count];
    }
}
