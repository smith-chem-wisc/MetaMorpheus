using Omics;
using Omics.Fragmentation;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class BioPolymerNotchFragmentIonComparer : Comparer<(int notch, IBioPolymerWithSetMods pwsm, List<MatchedFragmentIon> ions)>
    {
        /// <summary>
        /// Returns greater than 0 if x is better than y, less than 0 if y is better than x, and 0 if they are equal.
        /// Better is defined as having a lower notch, more fragment ions, or a shorter sequence (i.e. fewer modifications) in that order.
        /// If the aforementioned criteria are equal, then the two are compared based on the alphebetical ordering of the full sequence
        /// </summary>
        public override int Compare((int notch, IBioPolymerWithSetMods pwsm, List<MatchedFragmentIon> ions) x, (int notch, IBioPolymerWithSetMods pwsm, List<MatchedFragmentIon> ions) y)
        {
            if (x.notch != y.notch)
                return -1 * x.notch.CompareTo(y.notch); // Lower notch is better

            // Matched Ions is a nullable list, so we need to check for null
            if (x.ions == null && y.ions == null)
                return 0; // Both are null, they are equal
            if (x.ions == null)
                return -1; // x is null, y is better
            if (y.ions == null)
                return 1; // y is null, x is better
            if (x.ions.Count != y.ions.Count)
                return -1 * x.ions.Count.CompareTo(y.ions.Count); // More ions are better

            // Bpwsm is a nullable property, so we need to check for null
            if (x.pwsm == null && y.pwsm == null)
                return 0;
            if (x.pwsm == null)
                return 1; // Null Bpwsm is considered worse
            if (y.pwsm == null)
                return -1; // Null Bpwsm is considered worse

            if (x.pwsm.NumMods !=  y.pwsm.NumMods)
                return -1 * x.pwsm.NumMods.CompareTo(y.pwsm.NumMods); // Fewer mods are better

            if(x.pwsm.FullSequence != y.pwsm.FullSequence)
                return -1 * String.Compare(x.pwsm.FullSequence, y.pwsm.FullSequence); // (reverse) Alphabetical ordering of full sequence

            if(x.pwsm.Parent.Accession != y.pwsm.Parent.Accession) // This will break if the protein accession is not set (I'm not sure if that's possible)
                return -1 * String.Compare(x.pwsm.Parent.Accession, y.pwsm.Parent.Accession); // (reverse) Alphabetical ordering of protein accession

            return -1 * x.pwsm.OneBasedStartResidue.CompareTo(y.pwsm.OneBasedStartResidue);                                                                  
        }
    }
}
