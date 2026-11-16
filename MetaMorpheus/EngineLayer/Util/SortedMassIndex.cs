using System;

namespace EngineLayer.Util
{
    /// <summary>
    /// Finds the first entry at or above a value in an ascending double[] - the scan a precursor interval starts at, which
    /// ClassicSearchEngine looks up for every allowed interval of every peptide (66 intervals per peptide in the nightly
    /// GPTMD search) - by first reading where the value's fixed-width bucket starts, so that only a bucket's few entries are
    /// searched instead of the whole array. The answer is always the one
    /// <c>Array.BinarySearch</c> gives (its result, or its complement when the value is absent).
    /// </summary>
    internal sealed class SortedMassIndex
    {
        private readonly double[] Masses;

        /// <summary>BucketStarts[k] is the first index whose mass is at least k * BucketWidth. Null when the table is not used.</summary>
        private readonly int[] BucketStarts;

        /// <summary>A power of two, so that k * BucketWidth and value / BucketWidth are exact and the buckets meet exactly.</summary>
        private readonly double BucketWidth;

        /// <summary>The upper end of the last bucket; every mass in the array is below it.</summary>
        private readonly double TableEnd;

        public SortedMassIndex(double[] ascendingMasses)
        {
            Masses = ascendingMasses;

            // The table is only built over finite, non-negative masses in ascending order - what a sorted list of precursor
            // masses is. Anything else is left entirely to Array.BinarySearch.
            if (!IsAscendingFiniteAndNonNegative(ascendingMasses))
            {
                return;
            }

            // Half a dalton per bucket, widened (by doubling) until the table has at most about four buckets per mass, so
            // that a few wide-ranging masses cannot make it large.
            double maxMass = ascendingMasses[^1];
            double width = 0.5;
            int maxBuckets = Math.Max(1024, 4 * ascendingMasses.Length);
            while (maxMass / width >= maxBuckets)
            {
                width *= 2;
            }
            int bucketCount = (int)(maxMass / width) + 1;

            var bucketStarts = new int[bucketCount + 1];
            int j = 0;
            for (int k = 0; k <= bucketCount; k++)
            {
                double bucketStart = k * width;
                while (j < ascendingMasses.Length && ascendingMasses[j] < bucketStart)
                {
                    j++;
                }
                bucketStarts[k] = j;
            }

            BucketStarts = bucketStarts;
            BucketWidth = width;
            TableEnd = bucketCount * width;
        }

        /// <summary>
        /// The index of the first mass at or above <paramref name="value"/> (Masses.Length if there is none), exactly as
        /// <c>Array.BinarySearch</c> reports it: its result when it finds the value, otherwise the complement of its result.
        /// </summary>
        public int FirstAtOrAbove(double value)
        {
            // Beyond either end of the table no mass is equal, so BinarySearch's answer is its complement of the insertion point:
            // the length above (every mass is below TableEnd), 0 below (every mass is non-negative). Common, not an edge case:
            // the G-PTM-D search shifts every peptide's mass by each of its notches (66 in the bottom-up nightly's), and for a
            // heavy peptide many of the shifted intervals start above the heaviest precursor. Those fell through to
            // Array.BinarySearch over every scan.
            if (BucketStarts != null && value >= TableEnd)
            {
                return Masses.Length;
            }
            if (BucketStarts != null && value < 0)
            {
                return 0;
            }

            if (BucketStarts != null && value >= 0 && value < TableEnd)
            {
                // The answer lies between the start of the value's bucket and the start of the next one.
                int bucket = (int)(value / BucketWidth);
                int lo = BucketStarts[bucket];
                int hi = BucketStarts[bucket + 1];
                while (lo < hi)
                {
                    int mid = (int)((uint)(lo + hi) >> 1);
                    if (Masses[mid] < value)
                    {
                        lo = mid + 1;
                    }
                    else
                    {
                        hi = mid;
                    }
                }

                // An exact hit is left to BinarySearch: with repeated masses it may land on any of the equal entries.
                if (lo == Masses.Length || Masses[lo] != value)
                {
                    return lo;
                }
            }

            int index = Array.BinarySearch(Masses, value);
            return index < 0 ? ~index : index;
        }

        private static bool IsAscendingFiniteAndNonNegative(double[] masses)
        {
            if (masses.Length == 0)
            {
                return false;
            }
            for (int i = 0; i < masses.Length; i++)
            {
                // written so that NaN fails every test
                if (!(masses[i] >= 0) || double.IsPositiveInfinity(masses[i]) || (i > 0 && !(masses[i] >= masses[i - 1])))
                {
                    return false;
                }
            }
            return true;
        }
    }
}
