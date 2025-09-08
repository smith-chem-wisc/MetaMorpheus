using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Markup;
using MzLibUtil;

[assembly: System.Runtime.CompilerServices.InternalsVisibleTo("Test")]
namespace EngineLayer.Util
{
    /// <summary>
    /// Stores a set of unique <see cref="Precursor"/> instances, where uniqueness
    /// is defined by m/z within a specified tolerance.
    /// </summary>
    public class PrecursorSet : IEnumerable<Precursor>
    {
        private static ListPool<Precursor> _listPool = new(10);
        public readonly Tolerance Tolerance;
        /// <summary>
        /// This dictionary contains precursors indexed by their integer representation.
        /// The integer representation is calculated as the rounded m/z value multiplied by 100.
        /// </summary>
        internal Dictionary<int, List<Precursor>> PrecursorDictionary { get; }

        /// <summary>
        /// Initializes a new instance of the <see cref="PrecursorSet"/> class.
        /// </summary>
        /// <param name="mzTolerance">
        /// The m/z tolerance for determining if two precursors are considered equivalent.
        /// Must be positive.
        /// </param>
        public PrecursorSet(Tolerance tolerance)
        {
            Tolerance = tolerance;
            PrecursorDictionary = new Dictionary<int, List<Precursor>>();
        }

        public bool IsDirty { get; private set; } = false;

        /// <summary>
        /// Gets the number of precursors in the set.
        /// </summary>
        public int Count => PrecursorDictionary.Sum(kvp => kvp.Value.Count);

        /// <summary>
        /// Attempts to add a precursor to the set.
        /// </summary>
        /// <param name="precursor">The precursor to add.</param>
        /// <returns>
        /// True if the precursor was added; false if an equivalent precursor already exists.
        /// </returns>
        public bool Add(Precursor precursor)
        {
            if (precursor == null)
                return false;

            var integerKey = (int)Math.Round(precursor.MonoisotopicPeakMz * 100.0);
            if (!PrecursorDictionary.TryGetValue(integerKey, out var precursorsInBucket))
            {
                precursorsInBucket = new List<Precursor>();
                PrecursorDictionary[integerKey] = precursorsInBucket;
            }
            // Add the precursor to the appropriate bucket
            precursorsInBucket.Add(precursor);

            IsDirty = true;
            return true;
        }

        /// <summary>
        /// Checks whether the set contains an equivalent precursor within the specified m/z tolerance.
        /// </summary>
        /// <param name="precursor">The precursor to compare.</param>
        /// <returns>True if an equivalent precursor exists; otherwise, false.</returns>
        public bool ContainsEquivalent(Precursor precursor, out int integerKey)
        {
            if(precursor == null)
            {
                integerKey = 0;
                return false;
            }

            integerKey = (int)Math.Round(precursor.MonoisotopicPeakMz * 100.0);

            for (int i = integerKey - 1; i <= integerKey + 1; i++)
            {
                if(PrecursorDictionary.TryGetValue(i, out var precursorsInBucket))
                {
                    if (precursorsInBucket.Any(existingPrecursor => existingPrecursor.Equals(precursor, Tolerance)))
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        public void Sanitize()
        {
            if (!IsDirty)
                return;

            // Aggregate all precursors into a single list
            var allPrecursors = PrecursorDictionary.Values.SelectMany(list => list).ToList();

            MergeSplitEnvelopes(allPrecursors);
            FilterHarmonics(allPrecursors);
            RemoveDuplicates(allPrecursors); // Call this one last as it cleans and repopulates the dictionary

            IsDirty = false;
        }

        /// <summary>
        /// Classic decon sometimes splits a top-down peak across isotopic envelopes, this function merges them back together
        /// </summary>
        /// <param name="allPrecursors"></param>
        private void MergeSplitEnvelopes(List<Precursor> allPrecursors)
        {
            foreach (var chargeGroup in allPrecursors
                .Where(p => p.Envelope != null)
                .GroupBy(p => p.Charge))
            {
                
            }

        }

        private void FilterHarmonics(List<Precursor> allPrecursors)
        {

        }

        /// <summary>
        /// Currently only takes the first, future work would be to merge them
        /// </summary>
        private void RemoveDuplicates(List<Precursor> allPrecursors)
        {
            // Remove duplicates across all bins
            var uniquePrecursors = _listPool.Get();
            foreach (var precursor in allPrecursors)
            {
                bool isDuplicate = uniquePrecursors.Any(existing => precursor.Equals(existing, Tolerance));
                if (!isDuplicate)
                {
                    uniquePrecursors.Add(precursor);
                }
            }

            // Re-bin the unique precursors
            PrecursorDictionary.Clear();
            foreach (var precursor in uniquePrecursors)
            {
                var integerKey = (int)Math.Round(precursor.MonoisotopicPeakMz * 100.0);
                if (!PrecursorDictionary.TryGetValue(integerKey, out var bin))
                {
                    bin = new List<Precursor>();
                    PrecursorDictionary[integerKey] = bin;
                }
                bin.Add(precursor);
            }

            _listPool.Return(uniquePrecursors);
        }

        public void Clear()
        {
            foreach (var kvp in PrecursorDictionary)
            {
                kvp.Value.Clear();
            }
            PrecursorDictionary.Clear();
            IsDirty = false;
        }

        /// <summary>
        /// Returns an enumerator that iterates through the precursors in the set.
        /// </summary>
        public IEnumerator<Precursor> GetEnumerator()
        {
            if (IsDirty)
                Sanitize();
            return PrecursorDictionary.SelectMany(kvp => kvp.Value).GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
    }
}
