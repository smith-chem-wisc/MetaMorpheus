using Omics.Fragmentation;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    /// <summary>
    /// A peptide's theoretical product masses sorted once, ascending, with each one's index in the product list, so that
    /// every scan the peptide is scored against can find all of its closest envelopes in one forward walk of the scan's
    /// mass array (<see cref="Ms2ScanWithSpecificMass.GetClosestFragmentMasses"/>) instead of one binary search per product.
    /// Products without a mass (NaN) are left out; matching skips them anyway. Reused from peptide to peptide: call
    /// <see cref="Build"/> whenever the product list is regenerated.
    /// </summary>
    public sealed class SortedProductMasses
    {
        /// <summary>The non-NaN product masses, ascending; only the first <see cref="Count"/> are in use.</summary>
        internal double[] Masses = new double[64];

        /// <summary>For each entry of <see cref="Masses"/>, the index of its product in the list it was built from.</summary>
        internal int[] ProductIndices = new int[64];

        /// <summary>How many entries of <see cref="Masses"/> are in use.</summary>
        internal int Count;

        /// <summary>The length of the product list this was built from, NaN products included.</summary>
        internal int ProductCount;

        public void Build(List<Product> products)
        {
            if (Masses.Length < products.Count)
            {
                int size = Math.Max(products.Count, 2 * Masses.Length);
                Masses = new double[size];
                ProductIndices = new int[size];
            }

            int count = 0;
            for (int i = 0; i < products.Count; i++)
            {
                double mass = products[i].NeutralMass;
                if (!double.IsNaN(mass))
                {
                    Masses[count] = mass;
                    ProductIndices[count] = i;
                    count++;
                }
            }

            // Equal masses may come out in either order; that cannot matter, because the closest envelope of a product
            // depends on its mass alone.
            Array.Sort(Masses, ProductIndices, 0, count);
            Count = count;
            ProductCount = products.Count;
        }
    }
}
