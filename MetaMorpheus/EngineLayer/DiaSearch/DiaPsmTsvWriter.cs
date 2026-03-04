// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using MassSpectrometry.Dia;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace EngineLayer.DiaSearch
{
    /// <summary>
    /// Writes DIA search results to <c>.psmtsv</c> files that are compatible with MetaDraw
    /// and the standard MetaMorpheus PSM file format.
    ///
    /// Output format: 27 tab-separated columns.
    ///   Columns 1–11  match the MetaMorpheus .psmtsv standard exactly so MetaDraw's
    ///                 identification table populates correctly.
    ///   Columns 12–27 are DIA-specific extensions.
    ///
    /// File variants produced:
    ///   <c>AllDiaPSMs.psmtsv</c>                              — FDR-filtered aggregate across all files
    ///   <c>Individual File Results/*_DiaResults.psmtsv</c>    — FDR-filtered per file
    ///   <c>Individual File Results/*_PreFDR.psmtsv</c>        — all results (diagnostics)
    /// </summary>
    public static class DiaPsmTsvWriter
    {
        // ── Column headers ────────────────────────────────────────────────────

        /// <summary>
        /// Ordered column headers. First 11 are MetaMorpheus .psmtsv standard;
        /// columns 12-27 are DIA-specific.
        /// </summary>
        public static readonly IReadOnlyList<string> ColumnHeaders = new[]
        {
            // ── Standard MetaMorpheus .psmtsv columns (1-11) ─────────────────
            "File Name",
            "Scan Number",
            "Scan Retention Time",
            "Precursor Charge",
            "Precursor MZ",
            "Full Sequence",
            "Base Sequence",
            "Score",
            "QValue",
            "PEP",
            "PEP_QValue",

            // ── DIA-specific columns (12-27) ─────────────────────────────────
            "Peptide QValue",
            "Target/Decoy",
            "Fragments Detected",
            "Fragments Queried",
            "Fragment Detection Rate",
            "Spectral Angle Score",
            "Dot Product Score",
            "Apex Score",
            "Chimeric Score",
            "Classifier Score",
            "Library RT",
            "RT Window Start",
            "RT Window End",
            "RT Deviation (min)",
            "Modification Summary",
            "Window ID",
        };

        // ── Public API ────────────────────────────────────────────────────────

        /// <summary>
        /// Builds sorted adapter list from a raw results collection.
        /// Sorted by <see cref="DiaPsmAdapter.Score"/> descending,
        /// tie-broken by <see cref="DiaPsmAdapter.QValue"/> ascending.
        /// </summary>
        public static List<DiaPsmAdapter> BuildAdapters(
            IEnumerable<DiaSearchResult> results,
            string fileBaseName)
        {
            return results
                .Select(r => new DiaPsmAdapter(r, fileBaseName))
                .OrderByDescending(a => a.Score)
                .ThenBy(a => a.QValue)
                .ToList();
        }

        /// <summary>
        /// Writes a single-file .psmtsv. Rows are filtered by q-value threshold and
        /// optionally by target/decoy status.
        /// </summary>
        /// <param name="adapters">Adapters to write (pre-sorted by caller or BuildAdapters).</param>
        /// <param name="filePath">Full output path including filename and extension.</param>
        /// <param name="qValueThreshold">
        /// Maximum q-value to include. Pass <see cref="double.MaxValue"/> to write all rows.
        /// </param>
        /// <param name="includeDecoys">Whether to include decoy results.</param>
        public static void WriteToFile(
            IEnumerable<DiaPsmAdapter> adapters,
            string filePath,
            double qValueThreshold = 0.01,
            bool includeDecoys = false)
        {
            string dir = Path.GetDirectoryName(filePath);
            if (!string.IsNullOrEmpty(dir))
                Directory.CreateDirectory(dir);

            using var writer = new StreamWriter(filePath, false, Encoding.UTF8);
            WriteHeader(writer);

            foreach (var adapter in adapters)
            {
                if (!PassesFilter(adapter, qValueThreshold, includeDecoys))
                    continue;
                writer.WriteLine(FormatRow(adapter));
            }
        }

        /// <summary>
        /// Writes an aggregate .psmtsv from multiple per-file adapter collections.
        /// Rows across all files are merged, re-sorted by score, and filtered.
        /// </summary>
        /// <param name="perFileAdapters">One adapter list per input file.</param>
        /// <param name="filePath">Full output path.</param>
        /// <param name="qValueThreshold">Maximum q-value to include.</param>
        /// <param name="includeDecoys">Whether to include decoy results.</param>
        public static void WriteAggregate(
            IEnumerable<IEnumerable<DiaPsmAdapter>> perFileAdapters,
            string filePath,
            double qValueThreshold = 0.01,
            bool includeDecoys = false)
        {
            var merged = perFileAdapters
                .SelectMany(adapters => adapters)
                .OrderByDescending(a => a.Score)
                .ThenBy(a => a.QValue);

            string dir = Path.GetDirectoryName(filePath);
            if (!string.IsNullOrEmpty(dir))
                Directory.CreateDirectory(dir);

            using var writer = new StreamWriter(filePath, false, Encoding.UTF8);
            WriteHeader(writer);

            foreach (var adapter in merged)
            {
                if (!PassesFilter(adapter, qValueThreshold, includeDecoys))
                    continue;
                writer.WriteLine(FormatRow(adapter));
            }
        }

        // ── Internal helpers ──────────────────────────────────────────────────

        private static void WriteHeader(StreamWriter writer)
        {
            writer.WriteLine(string.Join("\t", ColumnHeaders));
        }

        private static bool PassesFilter(
            DiaPsmAdapter adapter,
            double qValueThreshold,
            bool includeDecoys)
        {
            if (adapter.UnderlyingResult.FdrInfo == null)
                return false;

            if (!includeDecoys && adapter.UnderlyingResult.IsDecoy)
                return false;

            if (adapter.QValue > qValueThreshold)
                return false;

            return true;
        }

        private static string FormatRow(DiaPsmAdapter a)
        {
            var r = a.UnderlyingResult;

            // Fragment detection rate (guard against div-by-zero)
            float fragDetRate = r.FragmentsQueried > 0
                ? (float)r.FragmentsDetected / r.FragmentsQueried
                : 0f;

            return string.Join("\t", new[]
            {
                // Standard columns 1-11
                a.FileName,
                a.ScanNumber.ToString(CultureInfo.InvariantCulture),
                a.ScanRetentionTime.ToString("F4", CultureInfo.InvariantCulture),
                a.PrecursorCharge.ToString(CultureInfo.InvariantCulture),
                a.PrecursorMz.ToString("F5", CultureInfo.InvariantCulture),
                a.FullSequence,
                a.BaseSequence,
                a.Score.ToString("F5", CultureInfo.InvariantCulture),
                a.QValue.ToString("E4", CultureInfo.InvariantCulture),
                FormatDouble(a.PEP),
                FormatDouble(a.PEP_QValue),

                // DIA-specific columns 12-27
                a.PeptideQValue.HasValue
                    ? a.PeptideQValue.Value.ToString("E4", CultureInfo.InvariantCulture)
                    : "",
                a.TargetDecoyLabel,
                a.FragmentsDetected.ToString(CultureInfo.InvariantCulture),
                a.FragmentsQueried.ToString(CultureInfo.InvariantCulture),
                fragDetRate.ToString("F4", CultureInfo.InvariantCulture),
                FormatFloat(a.SpectralAngleScore),
                FormatFloat(a.DotProductScore),
                FormatFloat(a.ApexScore),
                FormatFloat(a.ChimericScore),
                FormatFloat(r.ClassifierScore),
                r.LibraryRetentionTime.HasValue
                    ? r.LibraryRetentionTime.Value.ToString("F3", CultureInfo.InvariantCulture)
                    : "",
                a.RtWindowStart.ToString("F3", CultureInfo.InvariantCulture),
                a.RtWindowEnd.ToString("F3", CultureInfo.InvariantCulture),
                FormatFloat(a.RtDeviationMinutes),
                a.ModificationSummary,
                r.WindowId.ToString(CultureInfo.InvariantCulture),
            });
        }

        private static string FormatFloat(float value) =>
            float.IsNaN(value) ? "" : value.ToString("F5", CultureInfo.InvariantCulture);

        private static string FormatDouble(double value) =>
            double.IsNaN(value) ? "" : value.ToString("E4", CultureInfo.InvariantCulture);
    }
}