using System;

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
    /// was already amortizing. Measured on a synthetic version of the scoring loop, stamping runs
    /// ~8-13x faster than clearing when a scan touches ~0.1% of the index (a normal ppm-tolerance
    /// search) and ~0.6-0.8x as fast when it touches most of it (an open search). So the caller
    /// picks, and open searches keep clearing.
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
