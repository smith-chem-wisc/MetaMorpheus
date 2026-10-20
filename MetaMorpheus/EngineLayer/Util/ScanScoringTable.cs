using System;
using System.Linq;

namespace EngineLayer.Util
{
    /// <summary>
    /// The coarse (indexed) score of every peptide in the index, for the scan being scored.
    ///
    /// The search used to keep this as a byte[] and Array.Clear the whole thing between scans, which
    /// is O(peptides) of pure bookkeeping per scan whatever the scan actually touched. This can
    /// instead stamp each cell with the scan it belongs to and treat a stale stamp as a score of
    /// zero, which makes starting a scan O(1) -- the stamps only need clearing when they wrap, once
    /// every 255 scans.
    ///
    /// Stamping is not always the better trade. A stamped cell is two bytes wide, so a scan that
    /// touches most of the index pays double the memory traffic on the hot path to avoid a memset it
    /// was already amortizing. So the caller picks; see <see cref="IsWorthStamping"/>.
    ///
    /// Either way scores are bytes and increments wrap at 255, matching the byte[] this replaced.
    ///
    /// One instance per thread -- it is not safe to share.
    /// </summary>
    public sealed class ScanScoringTable
    {
        // Stamped layout: high byte is the scan stamp, low byte is the score.
        private readonly ushort[] _stampedCells;
        private readonly byte[] _scores;
        private byte _stamp;

        public ScanScoringTable(int peptideCount, bool stamped)
        {
            if (stamped)
            {
                _stampedCells = new ushort[peptideCount];
            }
            else
            {
                _scores = new byte[peptideCount];
            }
        }

        /// <summary>An unstamped table reading and writing <paramref name="scores"/> directly.</summary>
        internal ScanScoringTable(byte[] scores)
        {
            _scores = scores;
        }

        /// <summary>
        /// Whether a search with this acceptor should stamp rather than clear. Stamping pays off when a scan touches a
        /// small slice of the index, which is the case when the acceptor bounds the precursor mass on both sides: the
        /// coarse scoring then increments only the peptides inside that window. An acceptor unbounded on either side
        /// (OpenSearchMode, or ModOpen's [-187, +inf) interval) leaves the window open, so every bin is scored from
        /// its first peptide and a scan touches the index wholesale.
        /// </summary>
        public static bool IsWorthStamping(MassDiffAcceptor massDiffAcceptor)
        {
            // Whether a bound is infinite does not depend on the mass asked about, so any representative mass will do.
            const double representativePrecursorMass = 1000;
            return massDiffAcceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(representativePrecursorMass)
                .All(interval => !double.IsInfinity(interval.Minimum) && !double.IsInfinity(interval.Maximum));
        }

        /// <summary>
        /// The largest share of the peptide index a scan's precursor window may hold, on average, for stamping to be chosen.
        /// A bounded window is not necessarily a narrow one: a wide Custom interval is finite and can still hold most of
        /// the index, and it is how many cells a scan touches, not whether the window ends, that decides which table is faster.
        /// Measured end to end (modern search, 9,640 yeast MS2 scans against a 9.5M-peptide yeast index, 8 threads), stamping
        /// is faster up to a window holding about 2.5% of the index and breaks even near 5%; at 13% clearing is 1.14x faster.
        /// Ordinary tolerances hold under 0.1%. A synthetic scoring loop put the crossover at 0.2-0.3%, but its tables fit in
        /// cache; a real index's table does not, so each clear costs a full pass over memory.
        /// </summary>
        public const double MaxWindowShareWorthStamping = 0.05;

        /// <summary>
        /// Whether a search whose precursor windows hold, on average, <paramref name="meanWindowShareOfIndex"/> of the peptide
        /// index should stamp rather than clear. Only meaningful for an acceptor <see cref="IsWorthStamping"/> accepts.
        /// </summary>
        public static bool IsWindowNarrowEnoughToStamp(double meanWindowShareOfIndex)
        {
            return meanWindowShareOfIndex <= MaxWindowShareWorthStamping;
        }

        private bool Stamped => _stampedCells != null;

        /// <summary>Discards the previous scan's scores.</summary>
        public void BeginScan()
        {
            if (!Stamped)
            {
                Array.Clear(_scores, 0, _scores.Length);
                return;
            }

            if (_stamp == byte.MaxValue)
            {
                // Stamps have wrapped; retire every cell so stale ones cannot alias the new stamp.
                Array.Clear(_stampedCells, 0, _stampedCells.Length);
                _stamp = 0;
            }

            _stamp++;
        }

        public byte this[int peptideId]
        {
            get
            {
                if (!Stamped)
                {
                    return _scores[peptideId];
                }

                ushort cell = _stampedCells[peptideId];
                return (cell >> 8) == _stamp ? (byte)cell : (byte)0;
            }
        }

        public void Set(int peptideId, byte score)
        {
            if (!Stamped)
            {
                _scores[peptideId] = score;
                return;
            }

            _stampedCells[peptideId] = (ushort)((_stamp << 8) | score);
        }

        /// <summary>Adds one to a peptide's score and returns it, wrapping at 255 as a byte does.</summary>
        public byte Increment(int peptideId)
        {
            if (!Stamped)
            {
                return ++_scores[peptideId];
            }

            ushort cell = _stampedCells[peptideId];
            byte score = (cell >> 8) == _stamp ? (byte)((byte)cell + 1) : (byte)1;
            _stampedCells[peptideId] = (ushort)((_stamp << 8) | score);
            return score;
        }
    }
}
