#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;

namespace EngineLayer.ParallelSearch
{
    /// <summary>
    /// Per-file precursor and product mass tolerances for a search over scans merged from several files.
    /// Each scan carries the index of its file, so the single merged array and single engine pass are kept:
    /// the precursor window is searched at its widest and narrowed per scan, and the product tolerance is looked up per scan.
    /// </summary>
    public class PerFileTolerances
    {
        private PerFileTolerances(int[] scanFileIndex, Tolerance[] precursor, Tolerance[] product, Tolerance widestPrecursor)
        {
            ScanFileIndex = scanFileIndex;
            Precursor = precursor;
            Product = product;
            WidestPrecursor = widestPrecursor;
            PrecursorVaries = precursor.Any(p => !SameTolerance(p, precursor[0]));
            PrecursorIsPpm = widestPrecursor is PpmTolerance;

            if (PrecursorVaries)
            {
                // one read per candidate scan in the search: the half-width is centre * this (ppm) or this (absolute)
                double[] fileHalfWidth = precursor.Select(p => PrecursorIsPpm ? p.Value * 1e-6 : p.Value).ToArray();
                ScanPrecursorHalfWidth = new double[scanFileIndex.Length];
                for (int i = 0; i < scanFileIndex.Length; i++)
                {
                    ScanPrecursorHalfWidth[i] = fileHalfWidth[scanFileIndex[i]];
                }
            }
        }

        /// <summary>File index of each scan, aligned with the merged scan array.</summary>
        public int[] ScanFileIndex { get; }
        public Tolerance[] Precursor { get; }
        public Tolerance[] Product { get; }

        /// <summary>The widest precursor tolerance; build the acceptor from it so no file's candidates are missed.</summary>
        public Tolerance WidestPrecursor { get; }

        /// <summary>False when every file shares one precursor tolerance, so the per-scan filter can be skipped.</summary>
        public bool PrecursorVaries { get; }

        /// <summary>True when the precursor tolerances are ppm, so a scan's half-width scales with the theoretical mass.</summary>
        public bool PrecursorIsPpm { get; }

        /// <summary>
        /// Precursor half-width of each scan's file, aligned with the merged scan array: a relative factor for ppm,
        /// Daltons for absolute. Null when <see cref="PrecursorVaries"/> is false.
        /// </summary>
        public double[]? ScanPrecursorHalfWidth { get; }

        /// <summary>
        /// Returns null when every file uses the same tolerances as <paramref name="taskParams"/> (nothing to do),
        /// or when the tolerance types are mixed (ppm and absolute), which a single widest window cannot cover.
        /// </summary>
        public static PerFileTolerances? Create(IReadOnlyList<Ms2ScanWithSpecificMass> sortedScans,
            IEnumerable<(string FileName, CommonParameters Parameters)> fileParameters, CommonParameters taskParams,
            out string? warning)
        {
            warning = null;
            var byFile = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);
            var precursor = new List<Tolerance>();
            var product = new List<Tolerance>();
            foreach (var (fileName, parameters) in fileParameters)
            {
                if (byFile.ContainsKey(fileName)) continue;
                byFile[fileName] = precursor.Count;
                precursor.Add(parameters.PrecursorMassTolerance);
                product.Add(parameters.ProductMassTolerance);
            }

            if (precursor.All(p => SameTolerance(p, taskParams.PrecursorMassTolerance))
                && product.All(p => SameTolerance(p, taskParams.ProductMassTolerance)))
            {
                return null;
            }

            var fileIndex = new int[sortedScans.Count];
            for (int i = 0; i < fileIndex.Length; i++)
            {
                if (!byFile.TryGetValue(sortedScans[i].FullFilePath ?? string.Empty, out int index))
                {
                    // a scan from a file without an entry uses the task's tolerances
                    index = precursor.Count;
                    byFile[sortedScans[i].FullFilePath ?? string.Empty] = index;
                    precursor.Add(taskParams.PrecursorMassTolerance);
                    product.Add(taskParams.ProductMassTolerance);
                }
                fileIndex[i] = index;
            }

            // checked on the tolerances in use: the task's only counts when a scan's file has no entry of its own
            if (precursor.Select(p => p.GetType()).Distinct().Count() > 1)
            {
                warning = "File-specific precursor mass tolerances mix ppm and absolute units; using the task's precursor and product tolerances for every file.";
                return null;
            }

            Tolerance widest = precursor.OrderByDescending(p => p.Value).First();
            return new PerFileTolerances(fileIndex, precursor.ToArray(), product.ToArray(), widest);
        }

        private static bool SameTolerance(Tolerance a, Tolerance b) => a.GetType() == b.GetType() && a.Value == b.Value;
    }
}
