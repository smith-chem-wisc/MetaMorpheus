using System;

namespace EngineLayer
{
    /// <summary>
    /// Reusable working arrays for <see cref="MetaMorpheusEngine.ScoreFragmentIons"/>, one per thread, so that scoring a
    /// candidate peptide against a scan allocates nothing. Not thread-safe: give each thread its own.
    /// </summary>
    public sealed class FragmentMatchScratch
    {
        /// <summary>For each matched product, its index in the theoretical product list, in product order.</summary>
        internal int[] ProductIndices = new int[64];

        /// <summary>For each matched product, the index of its envelope in the scan's ExperimentalFragments.</summary>
        internal int[] EnvelopeIndices = new int[64];

        /// <summary>Closest envelope index of each sorted product mass, in sorted order.</summary>
        internal int[] SortedClosestIndices = new int[64];

        /// <summary>Closest envelope index of each product, by its index in the product list.</summary>
        internal int[] ClosestIndexByProduct = new int[64];

        internal void EnsureCapacity(int productCount)
        {
            if (ProductIndices.Length < productCount)
            {
                int size = Math.Max(productCount, 2 * ProductIndices.Length);
                ProductIndices = new int[size];
                EnvelopeIndices = new int[size];
                SortedClosestIndices = new int[size];
                ClosestIndexByProduct = new int[size];
            }
        }
    }
}
