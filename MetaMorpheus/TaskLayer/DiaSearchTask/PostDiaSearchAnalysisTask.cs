// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using EngineLayer;
using EngineLayer.DiaSearch;
using EngineLayer.DatabaseLoading;
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
            CalculateDiaFdr(Parameters.AllDiaResults, Parameters.ClassifierType);

            // ── Step 2: Filter and report ───────────────────────────────────────
            var passingResults = FilterResults(
                Parameters.AllDiaResults,
                Parameters.QValueThreshold,
                includeDecoys: Parameters.WriteDecoys,
                includeHighQValue: Parameters.WriteHighQValueResults);

            // ── Step 3: Write output files ──────────────────────────────────────
            // ── Step 3: Write output files ──────────────────────────────────────
            WriteDiaResults(passingResults, Parameters.OutputFolder, "AllDiaResults");

            if (Parameters.WriteIndividualFiles && Parameters.CurrentRawFileList?.Count > 1)
            {
                WritePerFileResults(Parameters.AllDiaResults, Parameters.CurrentRawFileList);
            }

            // ── Step 3b: Write .psmtsv output ───────────────────────────────────
            if (Parameters.WritePsmTsv)
            {
                WritePsmTsvResults(Parameters.AllDiaResults, Parameters.OutputFolder);
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
        /// Runs iterative semi-supervised FDR using the specified classifier
        /// (DiaFdrEngine.RunIterativeFdr from mzLib). The NeuralNetwork classifier
        /// produces 29,157 IDs at 1% FDR in the benchmark; LinearDiscriminant is
        /// ~4x faster with slightly fewer IDs and is used as fallback.
        /// </summary>
        private void CalculateDiaFdr(List<DiaSearchResult> results, DiaClassifierType classifierType)
        {
            Status("Computing feature vectors for iterative FDR...", Parameters.SearchTaskId);

            // ── Compute the 33-feature vectors ─────────────────────────────────
            DiaFeatureVector[] features = null;
            try
            {
                features = new DiaFeatureVector[results.Count];
                for (int i = 0; i < results.Count; i++)
                    features[i] = DiaFeatureExtractor.ComputeFeatures(results[i], i);
            }
            catch (Exception ex)
            {
                Warn($"Feature computation failed ({ex.Message}). Falling back to LDA FDR.");
            }

            // ── Iterative FDR (production path) ────────────────────────────────
            if (features != null && features.Length == results.Count)
            {
                int targetCount = results.Count(r => !r.IsDecoy);
                int decoyCount = results.Count(r => r.IsDecoy);

                Status($"Running iterative FDR on {results.Count:N0} results " +
                       $"({targetCount:N0} targets, {decoyCount:N0} decoys) with {classifierType}...",
                    Parameters.SearchTaskId);

                try
                {
                    var fdrResult = DiaFdrEngine.RunIterativeFdr(
                        results,
                        features,
                        classifierType: classifierType);

                    int idsAt1Pct = fdrResult.IdentificationsAt1PctFdr;
                    Status($"Iterative FDR complete: {idsAt1Pct:N0} IDs at 1% FDR " +
                           $"({fdrResult.IterationsCompleted} iterations, " +
                           $"classifier: {fdrResult.ClassifierType})",
                        Parameters.SearchTaskId);

                    return; // FdrInfo already set on each result by RunIterativeFdr
                }
                catch (Exception ex)
                {
                    Warn($"Iterative FDR with {classifierType} failed ({ex.Message}). " +
                         "Retrying with LinearDiscriminant fallback.");
                }

                // ── Fallback: retry with LDA if NN or GBT failed ───────────────
                if (classifierType != DiaClassifierType.LinearDiscriminant)
                {
                    try
                    {
                        Status("Running iterative FDR with LinearDiscriminant fallback...",
                            Parameters.SearchTaskId);
                        var fdrResult = DiaFdrEngine.RunIterativeFdr(
                            results,
                            features,
                            classifierType: DiaClassifierType.LinearDiscriminant);
                        Status($"LDA fallback FDR complete: {fdrResult.IdentificationsAt1PctFdr:N0} IDs at 1% FDR",
                            Parameters.SearchTaskId);
                        return;
                    }
                    catch (Exception ex2)
                    {
                        Warn($"LDA fallback FDR also failed ({ex2.Message}). Results will have no FDR annotation.");
                    }
                }
            }

            // ── Last resort: score-sort TDC without feature-based classifier ───
            Status("Running score-sort target-decoy FDR (last resort)...", Parameters.SearchTaskId);
            RunScoreSortFdr(results);
            Status("Score-sort FDR complete.", Parameters.SearchTaskId);
        }

        /// <summary>
        /// Minimal score-sort target-decoy competition. Used only when iterative FDR
        /// fails entirely. Sorts by SpectralAngleScore, walks TDC, assigns q-values.
        /// Does not populate PEP or peptide-level q-values.
        /// </summary>
        private static void RunScoreSortFdr(List<DiaSearchResult> results)
        {
            // Sort by spectral angle descending (best first)
            var sorted = results
                .Select((r, i) => (r, i))
                .OrderByDescending(x => x.r.SpectralAngleScore)
                .ToList();

            int cumTargets = 0, cumDecoys = 0;
            var qValues = new double[sorted.Count];

            for (int rank = 0; rank < sorted.Count; rank++)
            {
                if (sorted[rank].r.IsDecoy) cumDecoys++;
                else cumTargets++;
                qValues[rank] = cumTargets > 0 ? (double)cumDecoys / cumTargets : 1.0;
            }

            // Monotonize: walk back up, each q-value = min(q, next q)
            for (int rank = sorted.Count - 2; rank >= 0; rank--)
                qValues[rank] = Math.Min(qValues[rank], qValues[rank + 1]);

            for (int rank = 0; rank < sorted.Count; rank++)
            {
                var fdr = new DiaFdrInfo
                {
                    QValue = qValues[rank],
                    PEP = double.NaN,
                    PEP_QValue = double.NaN
                };
                sorted[rank].r.FdrInfo = fdr;
            }
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
        /// Writes DIA results in .psmtsv format for MetaDraw and downstream tools.
        /// Produces:
        ///   AllDiaPSMs.psmtsv                               — FDR-filtered aggregate
        ///   Individual File Results/*_DiaResults.psmtsv     — FDR-filtered per file
        ///   Individual File Results/*_DiaResults_PreFDR.psmtsv — all results (diagnostic)
        /// </summary>
        private void WritePsmTsvResults(List<DiaSearchResult> allResults, string outputFolder)
        {
            string fileBaseName = Parameters.CurrentRawFileList?.Count == 1
                ? Path.GetFileNameWithoutExtension(Parameters.CurrentRawFileList[0])
                : "DiaSearch";

            Status("Writing .psmtsv output...", Parameters.SearchTaskId);

            var adapters = DiaPsmTsvWriter.BuildAdapters(allResults, fileBaseName);

            // Per-file FDR-filtered .psmtsv
            string individualDir = Path.Combine(outputFolder, "Individual File Results");
            Directory.CreateDirectory(individualDir);

            DiaPsmTsvWriter.WriteToFile(
                adapters,
                Path.Combine(individualDir, $"{fileBaseName}_DiaResults.psmtsv"),
                qValueThreshold: Parameters.QValueThreshold,
                includeDecoys: Parameters.WriteDecoys);

            // Per-file pre-FDR diagnostic .psmtsv (all results, including decoys)
            DiaPsmTsvWriter.WriteToFile(
                adapters,
                Path.Combine(individualDir, $"{fileBaseName}_DiaResults_PreFDR.psmtsv"),
                qValueThreshold: double.MaxValue,
                includeDecoys: true);

            // Aggregate AllDiaPSMs.psmtsv
            DiaPsmTsvWriter.WriteAggregate(
                new[] { adapters },
                Path.Combine(outputFolder, "AllDiaPSMs.psmtsv"),
                qValueThreshold: Parameters.QValueThreshold,
                includeDecoys: Parameters.WriteDecoys);

            Status("Done writing .psmtsv output.", Parameters.SearchTaskId);
        }
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