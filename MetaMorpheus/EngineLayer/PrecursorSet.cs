using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;

[assembly: System.Runtime.CompilerServices.InternalsVisibleTo("Test")]
namespace EngineLayer
{
    /// <summary>
    /// Each precursor refers to one isotopic envelope observed in a single MS1 scan. 
    /// It is meant to hold experimentally observed data for a peptide/proteoform precursor.
    /// </summary>
    /// <param name="MonoisotopicPeakMz"></param>
    /// <param name="Charge"></param>
    /// <param name="Mass"></param>
    /// <param name="Intensity">Either the most abundant isotope intensity, or the summed intensity from all isotopes. For MM decon, this is determined by CommonParameters </param>
    /// <param name="EnvelopePeakCount"> The number of peaks observed in the Isotopic Envelope. 
    /// The minimum intensity for a peak to be considered is specified by the user in decon parameters, or in arguments passed to external software</param>
    /// <param name="FractionalIntensity"> The fraction of the intensity in the MS1 scan that is accounted for by this precursor</param>
    public record Precursor(double MonoisotopicPeakMz, int Charge, double Mass, double Intensity, int EnvelopePeakCount, double? FractionalIntensity = null);

    /// <summary>
    /// Stores a set of unique <see cref="Precursor"/> instances, where uniqueness
    /// is defined by m/z within a specified tolerance.
    /// </summary>
    public class PrecursorSet : IEnumerable<Precursor>
    {
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
            {
                return false;
            }

            if (ContainsEquivalent(precursor, out int integerKey))
            {
                return false;
            }

            if (!PrecursorDictionary.TryGetValue(integerKey, out var precursorsInBucket))
            {
                precursorsInBucket = new List<Precursor>();
                PrecursorDictionary[integerKey] = precursorsInBucket;
            }
            // Add the precursor to the appropriate bucket
            precursorsInBucket.Add(precursor);
            

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
                    foreach (var existingPrecursor in precursorsInBucket)
                    {
                        if (existingPrecursor.Charge == precursor.Charge
                            && Tolerance.Within(existingPrecursor.MonoisotopicPeakMz, precursor.MonoisotopicPeakMz))
                        {
                            return true;
                        }
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// Returns an enumerator that iterates through the precursors in the set.
        /// </summary>
        public IEnumerator<Precursor> GetEnumerator() => PrecursorDictionary.SelectMany(kvp => kvp.Value).GetEnumerator();

        IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();
    }
}
