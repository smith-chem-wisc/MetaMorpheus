using System;

namespace EngineLayer.Indexing
{
    /// <summary>
    /// The fragment index in compressed sparse row form: bin <c>b</c> holds the peptide ids
    /// <c>PeptideIds[BinStart[b] .. BinStart[b + 1])</c>.
    ///
    /// Replaces a <c>List&lt;int&gt;[]</c> of one small list object per occupied bin — several million of them
    /// for a mammalian database, each with its own header and repeatedly-doubled backing array. Two flat
    /// arrays remove that allocation and the garbage collection that went with it, and make the index
    /// writable and readable as raw bytes instead of object-by-object.
    ///
    /// Bins are empty rather than null, so callers test <c>IsEmpty</c> where they used to test for null. An
    /// empty bin and a null bin were already the same thing: bins were only ever created with an element in
    /// them.
    /// </summary>
    public sealed class FragmentIndex
    {
        private readonly int[] _binStart;
        private readonly int[] _peptideIds;

        /// <param name="binStart">Offsets into <paramref name="peptideIds"/>, length = bin count + 1.</param>
        /// <param name="peptideIds">Peptide ids, grouped by bin and ascending within each bin.</param>
        public FragmentIndex(int[] binStart, int[] peptideIds)
        {
            _binStart = binStart ?? throw new ArgumentNullException(nameof(binStart));
            _peptideIds = peptideIds ?? throw new ArgumentNullException(nameof(peptideIds));

            if (_binStart.Length < 1)
            {
                throw new ArgumentException("binStart needs at least one element", nameof(binStart));
            }
        }

        /// <summary>Number of bins.</summary>
        public int Length => _binStart.Length - 1;

        /// <summary>Total peptide ids across all bins.</summary>
        public int EntryCount => _peptideIds.Length;

        /// <summary>The peptide ids in a bin, ascending. Empty when nothing fragments into that bin.</summary>
        public ReadOnlySpan<int> this[int bin] => _peptideIds.AsSpan(_binStart[bin], _binStart[bin + 1] - _binStart[bin]);

        /// <summary>Bin offsets; exposed so the index can be written and read as raw arrays.</summary>
        public int[] BinStart => _binStart;

        /// <summary>The concatenated peptide ids.</summary>
        public int[] PeptideIds => _peptideIds;
    }
}
