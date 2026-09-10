using System;
using System.Collections.Generic;

namespace EngineLayer.Indexing
{
    /// <summary>
    /// The bins one contiguous block of peptides fragmented into, recorded in emission order.
    ///
    /// A block records its own emissions rather than writing into the shared index, which is what lets
    /// blocks fragment in parallel. Runs are concatenated in block order afterwards, so the sequence the
    /// builder sees is exactly the sequence a single sequential pass would have produced -- the index does
    /// not depend on how many blocks there were, or on which finished first.
    /// </summary>
    internal sealed class FragmentEmissionRun
    {
        private int[] _bins;
        private int[] _runLengths;
        private int _emissionCount;
        private int _peptideCount;
        private int _emissionsForCurrentPeptide;
        private bool _peptideOpen;

        internal FragmentEmissionRun(int firstPeptideId, int expectedPeptides = 16)
        {
            FirstPeptideId = firstPeptideId;
            _bins = new int[Math.Max(16, expectedPeptides)];
            _runLengths = new int[Math.Max(16, expectedPeptides)];
        }

        /// <summary>Peptide id the first run length belongs to. Ids ascend by one from here.</summary>
        internal int FirstPeptideId { get; }

        internal int[] Bins => _bins;
        internal int[] RunLengths => _runLengths;
        internal int EmissionCount => _emissionCount;
        internal int PeptideCount => _peptideCount;

        /// <summary>
        /// Starts the emissions for the next peptide id. Must be called once per peptide in the block, in
        /// ascending id order, including for peptides that emit nothing -- the run lengths are what map
        /// emissions back to peptide ids, so a skipped peptide would shift every id after it.
        /// </summary>
        internal void BeginPeptide()
        {
            if (_peptideOpen)
            {
                Append(ref _runLengths, _peptideCount, _emissionsForCurrentPeptide);
                _peptideCount++;
            }

            _peptideOpen = true;
            _emissionsForCurrentPeptide = 0;
        }

        internal void Add(int bin)
        {
            Append(ref _bins, _emissionCount, bin);
            _emissionCount++;
            _emissionsForCurrentPeptide++;
        }

        /// <summary>Closes the final peptide. Call once, after the last BeginPeptide/Add.</summary>
        internal void Complete()
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

    /// <summary>
    /// Turns per-block emission runs into a <see cref="FragmentIndex"/>.
    ///
    /// A stable counting sort over the concatenated runs: count how many emissions land in each bin,
    /// prefix-sum the counts into offsets, then scatter. Stability is what makes peptide ids ascend
    /// within a bin, which the search relies on -- BinarySearchBinForPrecursorIndex binary-searches a bin
    /// by peptide mass, and peptides are sorted by mass, so ascending id means ascending mass.
    /// </summary>
    internal static class FragmentIndexBuilder
    {
        internal static FragmentIndex Build(int binCount, IReadOnlyList<FragmentEmissionRun> runsInEmissionOrder)
        {
            // binStart doubles as the counting array and then as the scatter cursor, so a 30M-bin index
            // needs one 120 MB array here rather than three.
            int[] binStart;

            try
            {
                binStart = new int[binCount + 1];
            }
            catch (OutOfMemoryException)
            {
                throw new MetaMorpheusException("Max fragment mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
            }

            // Total in a long: the prefix sum below is int arithmetic and would wrap silently.
            long totalEntries = 0;
            foreach (FragmentEmissionRun run in runsInEmissionOrder)
            {
                int[] bins = run.Bins;
                int n = run.EmissionCount;
                for (int i = 0; i < n; i++)
                {
                    binStart[bins[i] + 1]++;
                }
                totalEntries += n;
            }

            if (totalEntries > int.MaxValue)
            {
                throw new MetaMorpheusException("Too many fragments to index; try \"Classic Search\" mode, more partitions, or a smaller maximum fragment mass");
            }

            for (int bin = 1; bin <= binCount; bin++)
            {
                binStart[bin] += binStart[bin - 1];
            }

            int[] peptideIds;
            try
            {
                peptideIds = new int[totalEntries];
            }
            catch (OutOfMemoryException)
            {
                throw new MetaMorpheusException("Max fragment mass too large for indexing engine; try \"Classic Search\" mode, or make the maximum fragment mass smaller");
            }

            // Scatter, using binStart itself as the per-bin write cursor.
            foreach (FragmentEmissionRun run in runsInEmissionOrder)
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
                        peptideIds[binStart[bins[emission++]]++] = peptideId;
                    }
                    peptideId++;
                }
            }

            // Scattering left binStart[bin] holding the start of bin + 1; shift it back.
            Array.Copy(binStart, 0, binStart, 1, binCount);
            binStart[0] = 0;

            return new FragmentIndex(binStart, peptideIds);
        }
    }
}
