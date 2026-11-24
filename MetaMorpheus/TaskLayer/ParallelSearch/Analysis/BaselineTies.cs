#nullable enable
using System.Collections.Generic;
using System.Linq;
using EngineLayer;

namespace TaskLayer.ParallelSearch.Analysis;

/// <summary>
/// Which transient matches tie with a protein outside the transient database, read before transient parsimony runs.
/// Parsimony trims each match's hypotheses in place, and when a transient protein explains a peptide that the only
/// baseline protein holding it also explains, it can keep the transient protein and drop the baseline one. The match
/// then lists only transient proteins and looks like the database's own evidence. On PXD061094 (ascites, sample 103),
/// HLA-A's FIAVGYVDDTQFVR lost P04439 this way in 22 PSMs of UP001521181, while the same PSMs kept it in sample 11.
/// </summary>
public sealed class BaselineTies
{
    private readonly HashSet<SpectralMatch> _tiedMatches;
    private readonly HashSet<string> _baselineSequences;

    private BaselineTies(HashSet<SpectralMatch> tiedMatches, HashSet<string> baselineSequences)
    {
        _tiedMatches = tiedMatches;
        _baselineSequences = baselineSequences;
    }

    public static BaselineTies Empty { get; } = new([], []);

    /// <summary>
    /// Reads the ties of these matches. A transient peptide whose sequence (I read as L) equals the sequence of a baseline
    /// hypothesis of the same match, target with target and decoy with decoy, is baseline-isobaric: it is remembered by
    /// sequence, so its other matches, and the peptide-level match, count as shared too.
    /// </summary>
    public static BaselineTies Read(IEnumerable<SpectralMatch> transientMatches, HashSet<string> transientProteinAccessions)
    {
        var tied = new HashSet<SpectralMatch>(ReferenceEqualityComparer.Instance);
        var sequences = new HashSet<string>();
        foreach (var match in transientMatches)
        {
            var hypotheses = match.HypothesesInAddedOrder.Select(h => h.SpecificBioPolymer).ToList();
            var baseline = hypotheses.Where(p => !transientProteinAccessions.Contains(p.Parent.Accession)).ToList();
            if (baseline.Count == 0)
                continue;
            tied.Add(match);
            foreach (var transient in hypotheses.Where(p => transientProteinAccessions.Contains(p.Parent.Accession)))
            {
                string key = Key(transient.BaseSequence, transient.Parent.IsDecoy);
                if (baseline.Any(b => Key(b.BaseSequence, b.Parent.IsDecoy) == key))
                    sequences.Add(key);
            }
        }
        return new BaselineTies(tied, sequences);
    }

    /// <summary>True when the match tied with a baseline protein before parsimony, or its transient peptide is baseline-isobaric.</summary>
    public bool Contains(SpectralMatch match, HashSet<string> transientProteinAccessions) =>
        _tiedMatches.Contains(match)
        || (_baselineSequences.Count > 0 && match.HypothesesInAddedOrder
            .Select(h => h.SpecificBioPolymer)
            .Any(p => transientProteinAccessions.Contains(p.Parent.Accession) && _baselineSequences.Contains(Key(p.BaseSequence, p.Parent.IsDecoy))));

    private static string Key(string baseSequence, bool isDecoy) => (isDecoy ? "D:" : "T:") + baseSequence.Replace('I', 'L');
}
