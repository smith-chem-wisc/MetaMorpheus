using System.Linq;
using EngineLayer.Util;

namespace EngineLayer.ClassicSearch
{
    /// <summary>
    /// The per-scan structures a classic search builds from its sorted scans: each scan's precursor mass for the search,
    /// the sorted index over those masses, and whether each scan can be scored before matching.
    /// <para>
    /// None of them depend on the proteins searched. ParallelSearch runs one classic search per transient database over
    /// the same scans, and rebuilt all of these for each one (about 60 bytes per scan per database, mostly on the large-object
    /// heap). Built once and passed to each <see cref="ClassicSearchEngine"/>, they are shared. Locks are not: searches
    /// running at once would wait on each other for the same scans, so each search keeps its own (striped) locks.
    /// </para>
    /// </summary>
    public sealed class ClassicSearchScanIndex
    {
        /// <param name="sortedScans">The scans every search using this index is given, the same array.</param>
        /// <param name="commonParameters">Parameters with the precursor-mass and fragment settings every such search uses.</param>
        public ClassicSearchScanIndex(Ms2ScanWithSpecificMass[] sortedScans, CommonParameters commonParameters)
        {
            Scans = sortedScans;
            PrecursorMasses = sortedScans.Select(b => b.GetPrecursorMassForSearch(commonParameters)).ToArray();
            MassIndex = new SortedMassIndex(PrecursorMasses);
            ScoreBeforeMatching = sortedScans
                .Select(scan => MetaMorpheusEngine.CanScoreBeforeMatching(scan, commonParameters, matchAllCharges: false))
                .ToArray();
        }

        public Ms2ScanWithSpecificMass[] Scans { get; }
        public double[] PrecursorMasses { get; }
        internal SortedMassIndex MassIndex { get; }

        /// <summary>Per scan, <see cref="MetaMorpheusEngine.CanScoreBeforeMatching"/> without matching all charges.</summary>
        public bool[] ScoreBeforeMatching { get; }

        /// <summary>Whether a search over <paramref name="scans"/> holding <paramref name="matchCount"/> matches can use this index.</summary>
        public bool Serves(Ms2ScanWithSpecificMass[] scans, int matchCount) =>
            ReferenceEquals(Scans, scans) && scans.Length == matchCount;

        /// <summary>The number of locks a search using this index stripes its scans over.</summary>
        public const int LockStripes = 4096;
    }
}
