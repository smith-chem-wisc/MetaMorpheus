using System;
using System.Buffers;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;

namespace EngineLayer.Indexing
{
    /// <summary>
    /// A fragment- or precursor-mass index, stored as one flat array of peptide ids plus a per-bin
    /// offset table (compressed sparse row).
    ///
    /// This replaces the List&lt;int&gt;[] the index used to be. That layout allocated a separate List
    /// object -- and its own backing int[], grown by doubling -- for every occupied bin, which meant
    /// millions of small objects for the GC to trace and a pointer chase on every bin lookup in the
    /// search inner loop. It also made the reference array alone cost 8 bytes per bin (240 MB at the
    /// default 30000 Da / 1000-bins-per-Da), which is what the OutOfMemoryException handler in
    /// IndexingEngine used to be apologizing for.
    ///
    /// Within a bin, peptide ids are ascending, exactly as appending in peptide order used to
    /// produce them. BinarySearchBinForPrecursorIndex relies on that: peptides are sorted by mass,
    /// so ascending peptide id within a bin means ascending mass within a bin.
    /// </summary>
    public sealed class MassBinIndex
    {
        // On-disk header. Bumping Version invalidates cached index files; the readers in
        // MetaMorpheusTask already fall back to regenerating when a read throws.
        private const int Magic = 0x4D4D4249; // "MMBI"
        private const int Version = 1;

        /// <summary>bin b occupies _entries[_binStarts[b].._binStarts[b + 1]). Length is BinCount + 1.</summary>
        private readonly int[] _binStarts;

        private readonly int[] _entries;

        private MassBinIndex(int[] binStarts, int[] entries)
        {
            _binStarts = binStarts;
            _entries = entries;
        }

        /// <summary>
        /// Number of mass bins. Named Length rather than Count because this stands in for the
        /// array it replaced, and callers bounds-check against it.
        /// </summary>
        public int Length => _binStarts.Length - 1;

        /// <summary>Total number of peptide ids across all bins.</summary>
        public int EntryCount => _entries.Length;

        /// <summary>The peptide ids in a bin, ascending. Empty -- not null -- when the bin is unoccupied.</summary>
        public ReadOnlySpan<int> this[int bin] => new ReadOnlySpan<int>(_entries, _binStarts[bin], _binStarts[bin + 1] - _binStarts[bin]);

        /// <summary>
        /// How many peptide ids are in a bin. Cheaper than materializing the span when the caller
        /// only wants to know whether the bin is worth visiting.
        /// </summary>
        public int CountInBin(int bin) => _binStarts[bin + 1] - _binStarts[bin];

        /// <summary>
        /// Builds the index from emission runs. Runs are concatenated in list order, so the caller
        /// controls the emission sequence; a stable counting sort over that sequence reproduces
        /// exactly the per-bin ordering that appending to List&lt;int&gt;[] used to produce.
        /// </summary>
        public static MassBinIndex Build(int binCount, IReadOnlyList<MassBinEmissionRun> runsInEmissionOrder)
        {
            // binStarts doubles as the counting array, then as the scatter cursor, so a 30M-bin
            // index needs one 120 MB array here rather than three.
            int[] binStarts = new int[binCount + 1];

            long total = 0;
            foreach (MassBinEmissionRun run in runsInEmissionOrder)
            {
                int[] bins = run.Bins;
                int n = run.EmissionCount;
                for (int i = 0; i < n; i++)
                {
                    binStarts[bins[i] + 1]++;
                }
                total += n;
            }

            if (total > int.MaxValue)
            {
                throw new MetaMorpheusException(
                    "Fragment index is too large to address (" + total + " entries); reduce the maximum fragment mass, " +
                    "raise the number of database partitions, or use Classic Search.");
            }

            for (int b = 0; b < binCount; b++)
            {
                binStarts[b + 1] += binStarts[b];
            }

            int[] entries = new int[(int)total];

            // Scatter using binStarts itself as the per-bin cursor. Afterwards binStarts[b] holds the
            // end of bin b, so shifting right by one turns the cursors back into starts.
            foreach (MassBinEmissionRun run in runsInEmissionOrder)
            {
                int[] bins = run.Bins;
                int[] runLengths = run.RunLengths;
                int peptideId = run.FirstPeptideId;
                int emission = 0;

                for (int p = 0; p < run.PeptideCount; p++)
                {
                    int emissionsForThisPeptide = runLengths[p];
                    for (int k = 0; k < emissionsForThisPeptide; k++)
                    {
                        entries[binStarts[bins[emission++]]++] = peptideId;
                    }
                    peptideId++;
                }
            }

            Array.Copy(binStarts, 0, binStarts, 1, binCount);
            binStarts[0] = 0;

            return new MassBinIndex(binStarts, entries);
        }

        public void Write(string filePath)
        {
            using FileStream file = File.Create(filePath);
            using BinaryWriter writer = new BinaryWriter(file);

            writer.Write(Magic);
            writer.Write(Version);
            writer.Write(_binStarts.Length);
            writer.Write(_entries.Length);

            WriteInts(file, _binStarts);
            WriteInts(file, _entries);
        }

        public static MassBinIndex Read(string filePath)
        {
            using FileStream file = File.OpenRead(filePath);
            using BinaryReader reader = new BinaryReader(file);

            if (reader.ReadInt32() != Magic || reader.ReadInt32() != Version)
            {
                // An index written by a different MetaMorpheus. Callers catch this and rebuild.
                throw new InvalidDataException("Index file is not in the format this version of MetaMorpheus writes.");
            }

            int binStartsLength = reader.ReadInt32();
            int entriesLength = reader.ReadInt32();

            if (binStartsLength < 1 || entriesLength < 0)
            {
                throw new InvalidDataException("Index file header is corrupt.");
            }

            int[] binStarts = new int[binStartsLength];
            int[] entries = new int[entriesLength];

            ReadInts(file, binStarts);
            ReadInts(file, entries);

            return new MassBinIndex(binStarts, entries);
        }

        // Copy through a rented buffer rather than one BinaryWriter.Write(int) per element; an index
        // can hold hundreds of millions of ints and the per-call overhead dominates otherwise.
        private const int CopyBufferBytes = 1 << 20;

        private static void WriteInts(Stream stream, int[] values)
        {
            byte[] buffer = ArrayPool<byte>.Shared.Rent(CopyBufferBytes);
            try
            {
                int intsPerChunk = buffer.Length / sizeof(int);
                for (int offset = 0; offset < values.Length; offset += intsPerChunk)
                {
                    int count = Math.Min(intsPerChunk, values.Length - offset);
                    Span<byte> bytes = buffer.AsSpan(0, count * sizeof(int));
                    MemoryMarshal.AsBytes(values.AsSpan(offset, count)).CopyTo(bytes);
                    stream.Write(bytes);
                }
            }
            finally
            {
                ArrayPool<byte>.Shared.Return(buffer);
            }
        }

        private static void ReadInts(Stream stream, int[] values)
        {
            byte[] buffer = ArrayPool<byte>.Shared.Rent(CopyBufferBytes);
            try
            {
                int intsPerChunk = buffer.Length / sizeof(int);
                for (int offset = 0; offset < values.Length; offset += intsPerChunk)
                {
                    int count = Math.Min(intsPerChunk, values.Length - offset);
                    Span<byte> bytes = buffer.AsSpan(0, count * sizeof(int));
                    stream.ReadExactly(bytes);
                    bytes.CopyTo(MemoryMarshal.AsBytes(values.AsSpan(offset, count)));
                }
            }
            finally
            {
                ArrayPool<byte>.Shared.Return(buffer);
            }
        }
    }

    /// <summary>
    /// An ordered run of (bin, peptide id) emissions covering a contiguous block of peptide ids.
    /// The peptide id is stored as a run length per peptide rather than once per emission, which
    /// halves the transient memory a build needs -- the fragment index of a large search can hold
    /// hundreds of millions of emissions.
    /// </summary>
    public sealed class MassBinEmissionRun
    {
        private int[] _bins;
        private int[] _runLengths;
        private int _emissionCount;
        private int _peptideCount;
        private int _emissionsForCurrentPeptide;
        private bool _peptideOpen;

        public MassBinEmissionRun(int firstPeptideId, int expectedPeptides = 16)
        {
            FirstPeptideId = firstPeptideId;
            _bins = new int[Math.Max(16, expectedPeptides)];
            _runLengths = new int[Math.Max(16, expectedPeptides)];
        }

        public int FirstPeptideId { get; }
        internal int[] Bins => _bins;
        internal int[] RunLengths => _runLengths;
        internal int EmissionCount => _emissionCount;
        internal int PeptideCount => _peptideCount;

        /// <summary>
        /// Starts the emissions for the next peptide id. Must be called once per peptide in the
        /// block, in ascending id order, including for peptides that emit nothing.
        /// </summary>
        public void BeginPeptide()
        {
            if (_peptideOpen)
            {
                Append(ref _runLengths, _peptideCount, _emissionsForCurrentPeptide);
                _peptideCount++;
            }

            _peptideOpen = true;
            _emissionsForCurrentPeptide = 0;
        }

        public void Add(int bin)
        {
            Append(ref _bins, _emissionCount, bin);
            _emissionCount++;
            _emissionsForCurrentPeptide++;
        }

        /// <summary>Closes the final peptide. Call once, after the last BeginPeptide/Add.</summary>
        public void Complete()
        {
            if (_peptideOpen)
            {
                Append(ref _runLengths, _peptideCount, _emissionsForCurrentPeptide);
                _peptideCount++;
                _peptideOpen = false;
            }
        }

        private static void Append(ref int[] array, int index, int value)
        {
            if (index == array.Length)
            {
                Array.Resize(ref array, array.Length * 2);
            }
            array[index] = value;
        }
    }
}
