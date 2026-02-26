// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using MassSpectrometry.Dia;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.FdrAnalysisDia
{
    /// <summary>
    /// False discovery rate analysis engine for DIA precursor matches.
    /// 
    /// Unlike the DDA FdrAnalysisEngine, this operates on DiaSearchResult objects which are
    /// fundamentally simpler: each result maps to exactly one library entry (target or decoy),
    /// has a single continuous score, and has no notch-level ambiguity or multi-protease grouping.
    /// 
    /// Algorithm:
    ///   1. Sort results by primary score (descending — best first)
    ///   2. Walk through sorted list counting cumulative targets and decoys
    ///   3. Compute raw FDR at each position: decoys / max(targets, 1)
    ///   4. Monotonize q-values by walking from worst to best (enforce non-increasing)
    ///   5. Repeat at peptide level (best score per unique sequence)
    /// 
    /// The engine modifies DiaSearchResult.FdrInfo in place, consistent with how the DDA
    /// FdrAnalysisEngine modifies SpectralMatch.FdrInfo.
    /// </summary>
    public class FdrAnalysisEngineDia : MetaMorpheusEngine
    {
        private readonly List<DiaSearchResult> _allResults;
        private readonly string _analysisType;
        private readonly DiaFdrScoreType _scoreType;

        /// <summary>
        /// Creates a new DIA FDR analysis engine.
        /// </summary>
        /// <param name="results">All DIA search results (targets + decoys). Will be sorted internally.</param>
        /// <param name="scoreType">Which score to use for ranking (DotProduct or SpectralAngle)</param>
        /// <param name="commonParameters">Common MetaMorpheus parameters (passed to base engine)</param>
        /// <param name="fileSpecificParameters">File-specific parameters (passed to base engine)</param>
        /// <param name="nestedIds">Task nesting IDs for status reporting</param>
        /// <param name="analysisType">Label for results reporting (default: "DIA Precursor")</param>
        public FdrAnalysisEngineDia(
            List<DiaSearchResult> results,
            DiaFdrScoreType scoreType,
            CommonParameters commonParameters,
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
            List<string> nestedIds,
            string analysisType = "DIA Precursor")
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            _allResults = results ?? throw new ArgumentNullException(nameof(results));
            _scoreType = scoreType;
            _analysisType = analysisType;

            // Ensure every result has a DiaFdrInfo object
            foreach (var result in _allResults)
            {
                result.FdrInfo ??= new DiaFdrInfo();
            }
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            var analysisResults = new FdrAnalysisResultsDia(this, _analysisType);

            Status("Running DIA FDR analysis...");

            if (GlobalVariables.StopLoops) return analysisResults;

            // Sort all results by score descending (best first)
            var sorted = SortByScore(_allResults, _scoreType);

            // ── Precursor-level q-values ────────────────────────────────────────
            CalculateQValues(sorted, peptideLevelCalculation: false);

            // ── Peptide-level q-values ──────────────────────────────────────────
            // Collapse to best score per unique sequence, recalculate, propagate back
            var peptides = CollapseToPeptides(sorted, _scoreType);
            CalculateQValues(peptides, peptideLevelCalculation: true);
            PropagatePeptideQValues(sorted, peptides);

            // ── Populate summary statistics ─────────────────────────────────────
            analysisResults.TotalTargets = sorted.Count(r => !r.IsDecoy);
            analysisResults.TotalDecoys = sorted.Count(r => r.IsDecoy);
            analysisResults.ResultsWithin1PercentFdr = sorted.Count(r =>
                !r.IsDecoy && r.FdrInfo.QValue <= 0.01);
            analysisResults.PeptidesWithin1PercentFdr = peptides.Count(r =>
                !r.IsDecoy && r.FdrInfo.PeptideQValue.HasValue && r.FdrInfo.PeptideQValue.Value <= 0.01);

            Status($"Done. {analysisResults.ResultsWithin1PercentFdr} precursors within 1% FDR, " +
                   $"{analysisResults.PeptidesWithin1PercentFdr} peptides within 1% FDR.");

            return analysisResults;
        }

        #region Core FDR Calculation

        /// <summary>
        /// Sorts results by the specified score type (descending).
        /// Results with NaN scores are placed last.
        /// </summary>
        internal static List<DiaSearchResult> SortByScore(List<DiaSearchResult> results, DiaFdrScoreType scoreType)
        {
            return results
                .OrderByDescending(r => GetScore(r, scoreType))
                .ThenByDescending(r => r.FragmentsDetected) // tiebreaker: more fragments detected
                .ToList();
        }

        /// <summary>
        /// Calculates q-values using standard target-decoy competition.
        /// 
        /// The list must already be sorted by score (descending — best first).
        /// 
        /// Algorithm:
        ///   1. Walk from best to worst score, counting cumulative targets and decoys.
        ///   2. Raw FDR at position i = cumDecoy / max(cumTarget, 1)
        ///   3. Monotonize: walk from worst to best, q[i] = min(q[i], q[i+1])
        ///      This ensures q-values are non-increasing as score improves,
        ///      matching the convention used in the DDA FdrAnalysisEngine.
        /// </summary>
        internal static void CalculateQValues(List<DiaSearchResult> sortedResults, bool peptideLevelCalculation)
        {
            if (sortedResults.Count == 0) return;

            double cumulativeTarget = 0;
            double cumulativeDecoy = 0;

            // ── Forward pass: compute raw FDR at each position ──────────────────
            for (int i = 0; i < sortedResults.Count; i++)
            {
                if (GlobalVariables.StopLoops) break;

                var result = sortedResults[i];

                if (result.IsDecoy)
                    cumulativeDecoy++;
                else
                    cumulativeTarget++;

                result.FdrInfo.CumulativeTarget = cumulativeTarget;
                result.FdrInfo.CumulativeDecoy = cumulativeDecoy;

                double rawFdr = cumulativeDecoy / Math.Max(cumulativeTarget, 1.0);
                double clampedFdr = Math.Min(rawFdr, 1.0);

                if (peptideLevelCalculation)
                    result.FdrInfo.PeptideQValue = clampedFdr;
                else
                    result.FdrInfo.QValue = clampedFdr;
            }

            // ── Backward pass: monotonize q-values ──────────────────────────────
            // Walk from the worst score (last) to the best (first), carrying the
            // running minimum. This ensures q-values never increase as score improves.
            double runningMin = 1.0;
            for (int i = sortedResults.Count - 1; i >= 0; i--)
            {
                if (peptideLevelCalculation)
                {
                    double current = sortedResults[i].FdrInfo.PeptideQValue!.Value;
                    runningMin = Math.Min(runningMin, current);
                    sortedResults[i].FdrInfo.PeptideQValue = runningMin;
                }
                else
                {
                    double current = sortedResults[i].FdrInfo.QValue;
                    runningMin = Math.Min(runningMin, current);
                    sortedResults[i].FdrInfo.QValue = runningMin;
                }
            }
        }

        /// <summary>
        /// Collapses results to peptide level: one result per unique sequence (best score wins).
        /// The returned list is sorted by score descending.
        /// Includes both targets and decoys for proper FDR calculation.
        /// </summary>
        internal static List<DiaSearchResult> CollapseToPeptides(
            List<DiaSearchResult> sortedResults, DiaFdrScoreType scoreType)
        {
            // Since sortedResults is already sorted by score descending,
            // taking the first per group gives the best-scoring result per sequence.
            return sortedResults
                .GroupBy(r => r.Sequence)
                .Select(g => g.First())
                .OrderByDescending(r => GetScore(r, scoreType))
                .ThenByDescending(r => r.FragmentsDetected)
                .ToList();
        }

        /// <summary>
        /// After peptide-level q-values are calculated on the collapsed list,
        /// propagates them back to all results sharing the same sequence.
        /// </summary>
        private static void PropagatePeptideQValues(
            List<DiaSearchResult> allResults, List<DiaSearchResult> peptideResults)
        {
            var peptideQValues = new Dictionary<string, double>(peptideResults.Count);
            foreach (var peptide in peptideResults)
            {
                if (peptide.FdrInfo.PeptideQValue.HasValue &&
                    !peptideQValues.ContainsKey(peptide.Sequence))
                {
                    peptideQValues[peptide.Sequence] = peptide.FdrInfo.PeptideQValue.Value;
                }
            }

            foreach (var result in allResults)
            {
                if (peptideQValues.TryGetValue(result.Sequence, out double pepQValue))
                {
                    result.FdrInfo.PeptideQValue = pepQValue;
                }
            }
        }

        #endregion

        #region Helpers

        /// <summary>
        /// Extracts the appropriate score from a DiaSearchResult.
        /// Returns float.NegativeInfinity for NaN scores so they sort to the bottom.
        /// </summary>
        private static float GetScore(DiaSearchResult result, DiaFdrScoreType scoreType)
        {
            float score = scoreType switch
            {
                DiaFdrScoreType.DotProduct => result.DotProductScore,
                DiaFdrScoreType.SpectralAngle => result.SpectralAngleScore,
                _ => result.SpectralAngleScore
            };

            return float.IsNaN(score) ? float.NegativeInfinity : score;
        }

        #endregion
    }
}
