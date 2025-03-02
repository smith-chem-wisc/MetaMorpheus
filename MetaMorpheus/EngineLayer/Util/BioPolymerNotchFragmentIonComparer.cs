using Omics;
using Omics.Fragmentation;
using System.Collections.Generic;

namespace EngineLayer;

/// <summary>
/// Compares the information of two tentative spectral matches to determine which is better
/// If used in an order by operation, the best matches will be first. 
/// </summary>
public class BioPolymerNotchFragmentIonComparer : Comparer<(int Notch, IBioPolymerWithSetMods Bpwsm, List<MatchedFragmentIon> MatchedIons)>
public class BioPolymerNotchFragmentIonComparer : Comparer<(int Notch, IBioPolymerWithSetMods Bpwsm, List<MatchedFragmentIon> MatchedIons)>,
    IComparer<TentativeSpectralMatch>
{
    /// <summary>
    /// Returns less than 0 if x is better than y, greater than 0 if y is better than x, and 0 if they are equal.
    /// Better is defined as having a lower Notch, more fragment ions, or a shorter sequence (i.e. fewer modifications) in that order.
    /// If the aforementioned criteria are equal, then the two are compared based on the alphebetical ordering of the full sequence
    /// </summary>
    public override int Compare((int Notch, IBioPolymerWithSetMods Bpwsm, List<MatchedFragmentIon> MatchedIons) x,
        (int Notch, IBioPolymerWithSetMods Bpwsm, List<MatchedFragmentIon> MatchedIons) y)
    {
        if (x.Notch != y.Notch)
            return x.Notch.CompareTo(y.Notch); // Lower Notch is better

        if (x.MatchedIons?.Count != y.MatchedIons?.Count && !ReferenceEquals(x.MatchedIons, null))
            return -1 * x.MatchedIons.Count.CompareTo(y.MatchedIons?.Count); // More ions are better

        if (x.Bpwsm == null && y.Bpwsm == null)
            return 0;
        if (x.Bpwsm == null)
            return 1; // Null Bpwsm is considered worse
        if (y.Bpwsm == null)
            return -1; // Null Bpwsm is considered worse

        if (x.Bpwsm.NumMods != y.Bpwsm.NumMods)
            return x.Bpwsm.NumMods.CompareTo(y.Bpwsm.NumMods); // Fewer mods are better

        if (x.Bpwsm.FullSequence != y.Bpwsm.FullSequence)
            return string.Compare(x.Bpwsm.FullSequence, y.Bpwsm.FullSequence); // Alphabetical ordering of full sequence

        if (x.Bpwsm.Parent?.Accession != y.Bpwsm.Parent?.Accession) // This will break if the protein accession is not set (I'm not sure if that's possible)
            return string.Compare(x.Bpwsm.Parent?.Accession, y.Bpwsm.Parent?.Accession); // Alphabetical ordering of protein accession

        return x.Bpwsm.OneBasedStartResidue.CompareTo(y.Bpwsm.OneBasedStartResidue);
    }

    public int Compare(TentativeSpectralMatch x, TentativeSpectralMatch y)
    {
        if (x is null && y is null)
            return 0;
        if (x is null) return 1;
        if (y is null) return -1;

        return Compare((x.Notch, x.WithSetMods, x.MatchedIons), (y.Notch, y.WithSetMods, y.MatchedIons));
    }
}