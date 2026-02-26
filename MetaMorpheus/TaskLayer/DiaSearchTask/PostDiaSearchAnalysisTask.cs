// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysisDia;
using MassSpectrometry.Dia;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace TaskLayer
{
    /// <summary>
    /// Post-search analysis task for DIA results.
    /// 
    /// Analogous to PostSearchAnalysisTask for DDA, but with a much simpler pipeline
    /// since DIA identification is library-based and doesn't require:
    /// - Protein parsimony (no protein inference from DIA XICs alone)
    /// - Modification localization (modifications come from the library)
    /// - FlashLFQ quantification (DIA quantification is XIC-based, already extracted)
    /// - Histogram analysis (no mass difference exploration in library search)
    /// 
    /// Pipeline:
    ///   1. Run FDR analysis (target-decoy competition, q-value calculation)
    ///   2. Filter results by q-value threshold
    ///   3. Write result files (all results TSV, per-file TSVs)
    ///   4. Report summary statistics
    /// </summary>
    public class PostDiaSearchAnalysisTask : MetaMorpheusTask
    {
        public PostDiaSearchAnalysisParameters Parameters { get; set; }

        public override string OutputFolder => Parameters?.OutputFolder;

        public PostDiaSearchAnalysisTask() : base(MyTask.Search)
        {
            // Uses MyTask.Search as base type (same as DiaSearchTask)
        }

        /// <summary>
        /// Main entry point — called by DiaSearchTask after extraction + scoring.
        /// </summary>
        public MyTaskResults Run()
        {
            if (GlobalVariables.StopLoops)
                return Parameters.SearchTaskResults;

            // ── Step 1: FDR Analysis ────────────────────────────────────────────
            CalculateDiaFdr(Parameters.AllDiaResults, Parameters.FdrScoreType);

            // ── Step 2: Filter and report ───────────────────────────────────────
            var passingResults = FilterResults(
                Parameters.AllDiaResults,
                Parameters.QValueThreshold,
                includeDecoys: Parameters.WriteDecoys,
                includeHighQValue: Parameters.WriteHighQValueResults);

            // ── Step 3: Write output files ──────────────────────────────────────
            WriteDiaResults(passingResults, Parameters.OutputFolder, "AllDiaResults");

            if (Parameters.WriteIndividualFiles && Parameters.CurrentRawFileList?.Count > 1)
            {
                WritePerFileResults(Parameters.AllDiaResults, Parameters.CurrentRawFileList);
            }

            // ── Step 4: Summary statistics ──────────────────────────────────────
            int targetsAt1Pct = Parameters.AllDiaResults.Count(r =>
                !r.IsDecoy && r.FdrInfo != null && r.FdrInfo.QValue <= 0.01);
            int peptidesAt1Pct = Parameters.AllDiaResults
                .Where(r => !r.IsDecoy && r.FdrInfo?.PeptideQValue != null && r.FdrInfo.PeptideQValue.Value <= 0.01)
                .Select(r => r.Sequence)
                .Distinct()
                .Count();

            Parameters.SearchTaskResults.AddTaskSummaryText(
                $"DIA precursors within 1% FDR: {targetsAt1Pct}");
            Parameters.SearchTaskResults.AddTaskSummaryText(
                $"DIA peptides within 1% FDR: {peptidesAt1Pct}");
            Parameters.SearchTaskResults.AddTaskSummaryText(
                $"Total target results: {Parameters.AllDiaResults.Count(r => !r.IsDecoy)}");
            Parameters.SearchTaskResults.AddTaskSummaryText(
                $"Total decoy results: {Parameters.AllDiaResults.Count(r => r.IsDecoy)}");

            return Parameters.SearchTaskResults;
        }

        /// <summary>
        /// Not used — DIA post-search is invoked via Run() directly.
        /// </summary>
        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            MyTaskResults = new MyTaskResults(this);
            return null;
        }

        #region FDR Calculation

        /// <summary>
        /// Runs FDR analysis on DIA results using the DIA-specific FDR engine.
        /// Modifies FdrInfo on each DiaSearchResult in place.
        /// </summary>
        private void CalculateDiaFdr(List<DiaSearchResult> results, DiaFdrScoreType scoreType)
        {
            Status("Estimating DIA FDR...", Parameters.SearchTaskId);

            var fdrEngine = new FdrAnalysisEngineDia(
                results,
                scoreType,
                CommonParameters,
                this.FileSpecificParameters,
                new List<string> { Parameters.SearchTaskId },
                analysisType: "DIA Precursor");

            fdrEngine.Run();

            Status("Done estimating DIA FDR!", Parameters.SearchTaskId);
        }

        #endregion

        #region Result Filtering

        /// <summary>
        /// Filters DIA results based on q-value threshold and inclusion settings.
        /// Returns a new list (does not modify the original).
        /// </summary>
        internal static List<DiaSearchResult> FilterResults(
            List<DiaSearchResult> allResults,
            double qValueThreshold,
            bool includeDecoys = false,
            bool includeHighQValue = false)
        {
            return allResults.Where(r =>
            {
                // Skip results without FDR info
                if (r.FdrInfo == null) return false;

                // Filter decoys unless requested
                if (r.IsDecoy && !includeDecoys) return false;

                // Filter by q-value unless writing all
                if (!includeHighQValue && r.FdrInfo.QValue > qValueThreshold) return false;

                return true;
            }).ToList();
        }

        #endregion

        #region Result Output

        /// <summary>
        /// Writes DIA results to a tab-separated file.
        /// </summary>
        private void WriteDiaResults(List<DiaSearchResult> results, string outputFolder, string fileLabel)
        {
            string filePath = Path.Combine(outputFolder, fileLabel + ".tsv");

            Status($"Writing DIA results to {fileLabel}.tsv...", Parameters.SearchTaskId);

            using var writer = new StreamWriter(filePath);

            // Header
            writer.WriteLine(GetTsvHeader());

            // Data rows — sorted by q-value ascending (best first)
            var sorted = results
                .OrderBy(r => r.FdrInfo?.QValue ?? 2.0)
                .ThenByDescending(r => r.SpectralAngleScore);

            foreach (var result in sorted)
            {
                writer.WriteLine(GetTsvLine(result));
            }

            Status($"Done writing {results.Count} results to {fileLabel}.tsv", Parameters.SearchTaskId);
        }

        /// <summary>
        /// Writes per-file DIA results when multiple raw files were searched.
        /// </summary>
        private void WritePerFileResults(List<DiaSearchResult> allResults, List<string> rawFileList)
        {
            // For now, DIA results don't carry a source file reference
            // (all precursors are from the library, matched against whichever file was searched).
            // Per-file output will be implemented when multi-file DIA search is supported.
            // This is a placeholder for the architecture.
        }

        /// <summary>
        /// Returns the TSV header line for DIA results output.
        /// </summary>
        internal static string GetTsvHeader()
        {
            return string.Join("\t", new[]
            {
                "Sequence",
                "Charge",
                "Precursor m/z",
                "Window ID",
                "Is Decoy",
                "Dot Product Score",
                "Spectral Angle Score",
                "Fragments Detected",
                "Fragments Queried",
                "Fragment Detection Rate",
                "Library RT",
                "RT Window Start",
                "RT Window End",
                "QValue",
                "Peptide QValue",
                "PEP",
                "PEP QValue"
            });
        }

        /// <summary>
        /// Converts a single DIA result to a TSV data line.
        /// </summary>
        internal static string GetTsvLine(DiaSearchResult result)
        {
            var fdr = result.FdrInfo;

            return string.Join("\t", new[]
            {
                result.Sequence,
                result.ChargeState.ToString(CultureInfo.InvariantCulture),
                result.PrecursorMz.ToString("F5", CultureInfo.InvariantCulture),
                result.WindowId.ToString(CultureInfo.InvariantCulture),
                result.IsDecoy ? "TRUE" : "FALSE",
                FormatFloat(result.DotProductScore),
                FormatFloat(result.SpectralAngleScore),
                result.FragmentsDetected.ToString(CultureInfo.InvariantCulture),
                result.FragmentsQueried.ToString(CultureInfo.InvariantCulture),
                result.FragmentDetectionRate.ToString("F4", CultureInfo.InvariantCulture),
                result.LibraryRetentionTime?.ToString("F3", CultureInfo.InvariantCulture) ?? "",
                result.RtWindowStart.ToString("F3", CultureInfo.InvariantCulture),
                result.RtWindowEnd.ToString("F3", CultureInfo.InvariantCulture),
                fdr != null ? fdr.QValue.ToString("E4", CultureInfo.InvariantCulture) : "",
                fdr?.PeptideQValue?.ToString("E4", CultureInfo.InvariantCulture) ?? "",
                fdr != null && !double.IsNaN(fdr.PEP) ? fdr.PEP.ToString("E4", CultureInfo.InvariantCulture) : "",
                fdr != null && !double.IsNaN(fdr.PEP_QValue) ? fdr.PEP_QValue.ToString("E4", CultureInfo.InvariantCulture) : ""
            });
        }

        private static string FormatFloat(float value)
        {
            return float.IsNaN(value) ? "" : value.ToString("F5", CultureInfo.InvariantCulture);
        }

        #endregion
    }
}