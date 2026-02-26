// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;

namespace Development.Dia
{
    /// <summary>
    /// End-to-end benchmark and validation for the full DIA pipeline:
    ///   Synthetic DIA scans → Index → Query Gen → Extraction → Scoring → FDR → TSV
    ///
    /// Generates entirely synthetic data — no external files needed.
    ///
    /// Design:
    ///   - Target library entries have fragments planted into scan data at known m/z + RT
    ///   - Decoy library entries use +50 Da shifted fragments that miss planted signals
    ///   - Targets should score high, decoys low → FDR correctly separates them
    ///
    /// Run via Development project: Development.Program.Main → DiaEndToEndBenchmark.RunAll()
    /// </summary>
    public static class DiaEndToEndBenchmark
    {
        public static void RunAll()
        {
            Console.WriteLine("=== DIA End-to-End Pipeline Benchmark ===");
            Console.WriteLine();

            RunSmallValidation();
            Console.WriteLine(new string('─', 60));
            RunLargeScale();
            Console.WriteLine(new string('─', 60));
            RunEdgeCases();
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Benchmark 1: Small-scale validation (correctness focus)
        // ──────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Small dataset focused on correctness validation.
        /// 10 windows, 100 scans/window, 50 targets, 50 decoys.
        /// Verifies: scores, q-value monotonicity, T/D separation, TSV output.
        /// </summary>
        public static void RunSmallValidation()
        {
            Console.WriteLine("── Small Validation (correctness) ──");
            var totalSw = Stopwatch.StartNew();

            const int windowCount = 10;
            const int scansPerWindow = 100;
            const int peaksPerScan = 200;
            const int targetCount = 50;
            const int decoyCount = 50;
            const int fragmentsPerPrecursor = 6;
            var rng = new Random(42);

            // Generate
            var sw = Stopwatch.StartNew();
            var (scans, targets, decoys) = GenerateSyntheticDiaData(
                windowCount, scansPerWindow, peaksPerScan,
                targetCount, decoyCount, fragmentsPerPrecursor,
                400.0, 20.0, rng);
            sw.Stop();
            Console.WriteLine($"  Data gen:        {sw.ElapsedMilliseconds,6} ms  ({scans.Length} scans)");

            // Build index
            sw.Restart();
            using var index = DiaScanIndexBuilder.Build(scans);
            sw.Stop();
            Console.WriteLine($"  Index build:     {sw.ElapsedMilliseconds,6} ms  ({index.WindowCount} windows, {index.ScanCount} scans)");

            // Generate queries
            var allPrecursors = new List<LibraryPrecursorInput>(targets.Count + decoys.Count);
            allPrecursors.AddRange(targets);
            allPrecursors.AddRange(decoys);

            var parameters = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 5.0f,
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0.0f,
                MaxThreads = 4,
                PreferGpu = false
            };

            sw.Restart();
            var genResult = DiaLibraryQueryGenerator.Generate(allPrecursors, index, parameters);
            sw.Stop();
            Console.WriteLine($"  Query gen:       {sw.ElapsedMilliseconds,6} ms  ({genResult.Queries.Length:N0} queries, " +
                $"{genResult.PrecursorGroups.Length} precursors, " +
                $"{genResult.SkippedNoWindow} skipped-window, {genResult.SkippedNoFragments} skipped-frag)");

            // Extract
            sw.Restart();
            using var orchestrator = new DiaExtractionOrchestrator(index);
            var extraction = orchestrator.ExtractAll(genResult.Queries, maxDegreeOfParallelism: 4);
            sw.Stop();
            Console.WriteLine($"  Extraction:      {sw.ElapsedMilliseconds,6} ms  ({extraction.TotalDataPoints:N0} data points)");

            // Score + assemble
            sw.Restart();
            var dotScorer = new NormalizedDotProductScorer();
            var saScorer = new SpectralAngleScorer();
            var diaResults = DiaLibraryQueryGenerator.AssembleResults(
                allPrecursors, genResult, extraction.Results, parameters,
                dotProductScorer: dotScorer, spectralAngleScorer: saScorer);
            sw.Stop();

            var targetResults = diaResults.Where(r => !r.IsDecoy).ToList();
            var decoyResults = diaResults.Where(r => r.IsDecoy).ToList();
            Console.WriteLine($"  Score+assemble:  {sw.ElapsedMilliseconds,6} ms  ({diaResults.Count} results: " +
                $"{targetResults.Count} T, {decoyResults.Count} D)");

            // FDR
            sw.Restart();
            var fdrResults = RunFdr(diaResults);
            sw.Stop();
            Console.WriteLine($"  FDR:             {sw.ElapsedMilliseconds,6} ms");

            // TSV output
            string outputDir = Path.Combine(Path.GetTempPath(), "DiaE2E_" + Guid.NewGuid().ToString("N")[..8]);
            Directory.CreateDirectory(outputDir);
            sw.Restart();
            string tsvPath = WriteTsv(diaResults, Path.Combine(outputDir, "AllDiaResults.tsv"));
            sw.Stop();
            int tsvLines = File.ReadAllLines(tsvPath).Length;
            Console.WriteLine($"  TSV write:       {sw.ElapsedMilliseconds,6} ms  ({tsvLines - 1} rows → {tsvPath})");

            totalSw.Stop();
            Console.WriteLine($"  TOTAL:           {totalSw.ElapsedMilliseconds,6} ms");
            Console.WriteLine();

            // ── Validation checks ──────────────────────────────────────────
            int passed = 0;
            int failed = 0;

            // Check 1: Targets score higher on average
            if (targetResults.Count > 0 && decoyResults.Count > 0)
            {
                double avgTargetSA = targetResults.Where(r => !float.IsNaN(r.SpectralAngleScore))
                    .DefaultIfEmpty().Average(r => r?.SpectralAngleScore ?? 0);
                double avgDecoySA = decoyResults.Where(r => !float.IsNaN(r.SpectralAngleScore))
                    .DefaultIfEmpty().Average(r => r?.SpectralAngleScore ?? 0);

                Check(ref passed, ref failed,
                    avgTargetSA > avgDecoySA,
                    $"Target avg SA ({avgTargetSA:F4}) > Decoy avg SA ({avgDecoySA:F4})");
            }

            // Check 2: Q-values monotonically non-decreasing (best→worst score)
            var sorted = diaResults
                .OrderByDescending(r => float.IsNaN(r.SpectralAngleScore) ? -1f : r.SpectralAngleScore)
                .ToList();
            bool monotonic = true;
            for (int i = 1; i < sorted.Count; i++)
            {
                if (sorted[i].FdrInfo.QValue < sorted[i - 1].FdrInfo.QValue - 1e-10)
                {
                    monotonic = false;
                    break;
                }
            }
            Check(ref passed, ref failed, monotonic, "Q-values monotonically non-decreasing");

            // Check 3: Some targets pass 1% FDR
            int at1Pct = targetResults.Count(r => r.FdrInfo.QValue <= 0.01);
            Check(ref passed, ref failed, at1Pct > 0, $"Targets at 1% FDR: {at1Pct}");

            // Check 4: All results have FdrInfo
            bool allHaveFdr = diaResults.All(r => r.FdrInfo != null);
            Check(ref passed, ref failed, allHaveFdr, "All results have FdrInfo");

            // Check 5: All results have peptide q-values
            bool allHavePepQv = diaResults.All(r => r.FdrInfo?.PeptideQValue != null);
            Check(ref passed, ref failed, allHavePepQv, "All results have peptide q-values");

            // Check 6: Q-values in valid range [0, 1]
            bool qvRange = diaResults.All(r => r.FdrInfo.QValue >= 0 && r.FdrInfo.QValue <= 1.0);
            Check(ref passed, ref failed, qvRange, "All q-values in [0, 1]");

            // Check 7: TSV header has expected columns
            string header = File.ReadLines(tsvPath).First();
            bool hasQv = header.Contains("QValue");
            bool hasPepQv = header.Contains("Peptide QValue");
            bool hasSA = header.Contains("Spectral Angle Score");
            Check(ref passed, ref failed, hasQv && hasPepQv && hasSA,
                "TSV has QValue, Peptide QValue, Spectral Angle Score columns");

            // Check 8: TSV data rows parseable
            var dataLines = File.ReadAllLines(tsvPath);
            bool parseable = true;
            if (dataLines.Length > 1)
            {
                string[] hdCols = dataLines[0].Split('\t');
                int qvIdx = Array.IndexOf(hdCols, "QValue");
                if (qvIdx >= 0)
                {
                    string[] row = dataLines[1].Split('\t');
                    parseable = row.Length == hdCols.Length &&
                        double.TryParse(row[qvIdx], NumberStyles.Float, CultureInfo.InvariantCulture, out _);
                }
                else parseable = false;
            }
            Check(ref passed, ref failed, parseable, "TSV data rows parseable");

            Console.WriteLine($"  Validation: {passed} passed, {failed} failed");
            Console.WriteLine();

            // Cleanup
            try { Directory.Delete(outputDir, true); } catch { }
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Benchmark 2: Large-scale performance
        // ──────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Larger dataset focused on throughput measurement.
        /// 20 windows, 500 scans/window, 500 targets, 500 decoys, 8 fragments each.
        /// </summary>
        public static void RunLargeScale()
        {
            Console.WriteLine("── Large-Scale Benchmark (throughput) ──");
            var rng = new Random(42);

            const int windowCount = 20;
            const int scansPerWindow = 500;
            const int peaksPerScan = 300;
            const int targetCount = 500;
            const int decoyCount = 500;
            const int fragmentsPerPrecursor = 8;

            long totalPeaks = (long)windowCount * scansPerWindow * peaksPerScan;
            Console.WriteLine($"  Config: {windowCount} win × {scansPerWindow} scans × {peaksPerScan} peaks = " +
                $"{totalPeaks:N0} peaks");
            Console.WriteLine($"          {targetCount} targets + {decoyCount} decoys × {fragmentsPerPrecursor} frags");
            Console.WriteLine();

            var totalSw = Stopwatch.StartNew();

            var sw = Stopwatch.StartNew();
            var (scans, targets, decoys) = GenerateSyntheticDiaData(
                windowCount, scansPerWindow, peaksPerScan,
                targetCount, decoyCount, fragmentsPerPrecursor,
                400.0, 20.0, rng);
            sw.Stop();
            Console.WriteLine($"  Data gen:        {sw.ElapsedMilliseconds,6} ms");

            sw.Restart();
            using var index = DiaScanIndexBuilder.Build(scans);
            sw.Stop();
            Console.WriteLine($"  Index build:     {sw.ElapsedMilliseconds,6} ms  ({totalPeaks / Math.Max(1, sw.Elapsed.TotalSeconds):N0} peaks/sec)");

            var allPrecursors = targets.Concat(decoys).ToList();
            var parameters = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 5.0f,
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0.0f,
                MaxThreads = -1, // use all cores
                PreferGpu = false
            };

            sw.Restart();
            var genResult = DiaLibraryQueryGenerator.Generate(allPrecursors, index, parameters);
            sw.Stop();
            Console.WriteLine($"  Query gen:       {sw.ElapsedMilliseconds,6} ms  ({genResult.Queries.Length:N0} queries)");

            sw.Restart();
            using var orchestrator = new DiaExtractionOrchestrator(index);
            var extraction = orchestrator.ExtractAll(genResult.Queries);
            sw.Stop();
            double queriesPerSec = genResult.Queries.Length / Math.Max(0.001, sw.Elapsed.TotalSeconds);
            Console.WriteLine($"  Extraction:      {sw.ElapsedMilliseconds,6} ms  ({queriesPerSec:N0} queries/sec)");

            sw.Restart();
            var diaResults = DiaLibraryQueryGenerator.AssembleResults(
                allPrecursors, genResult, extraction.Results, parameters,
                new NormalizedDotProductScorer(), new SpectralAngleScorer());
            sw.Stop();
            Console.WriteLine($"  Score+assemble:  {sw.ElapsedMilliseconds,6} ms  ({diaResults.Count} results)");

            sw.Restart();
            var fdrResults = RunFdr(diaResults);
            sw.Stop();
            int at1Pct = diaResults.Count(r => !r.IsDecoy && r.FdrInfo.QValue <= 0.01);
            Console.WriteLine($"  FDR:             {sw.ElapsedMilliseconds,6} ms  ({at1Pct} targets at 1%)");

            totalSw.Stop();
            Console.WriteLine();
            Console.WriteLine($"  TOTAL:           {totalSw.ElapsedMilliseconds,6} ms");
            Console.WriteLine($"  Processors:      {Environment.ProcessorCount}");
            Console.WriteLine();
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Edge case validation
        // ──────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Tests FDR edge cases: no decoys, equal scores, peptide-level collapsing.
        /// </summary>
        public static void RunEdgeCases()
        {
            Console.WriteLine("── Edge Case Validation ──");
            int passed = 0;
            int failed = 0;

            // Edge 1: No decoys → all q-values = 0
            {
                var results = new List<DiaSearchResult>();
                for (int i = 0; i < 20; i++)
                {
                    var r = new DiaSearchResult($"PEPTIDE{i}", 2, 500.0 + i, 0, false, 6, 30.0, 25f, 35f);
                    r.SpectralAngleScore = 0.9f - i * 0.01f;
                    r.DotProductScore = 0.85f - i * 0.01f;
                    r.FragmentsDetected = 5;
                    results.Add(r);
                }
                RunFdr(results);

                bool allZero = results.All(r => r.FdrInfo.QValue == 0.0);
                Check(ref passed, ref failed, allZero, "No decoys → all q-values = 0");
            }

            // Edge 2: Equal T/D scores → high FDR, few pass 1%
            {
                var results = new List<DiaSearchResult>();
                var rng = new Random(123);
                for (int i = 0; i < 50; i++)
                {
                    float score = (float)rng.NextDouble();

                    var t = new DiaSearchResult($"TARGET{i}", 2, 500.0 + i, 0, false, 6, 30.0, 25f, 35f);
                    t.SpectralAngleScore = score;
                    t.FragmentsDetected = 5;
                    results.Add(t);

                    var d = new DiaSearchResult($"DECOY{i}", 2, 500.0 + i, 0, true, 6, 30.0, 25f, 35f);
                    d.SpectralAngleScore = score + (float)(rng.NextDouble() * 0.02 - 0.01);
                    d.FragmentsDetected = 5;
                    results.Add(d);
                }
                RunFdr(results);

                int atOnePct = results.Count(r => !r.IsDecoy && r.FdrInfo.QValue <= 0.01);
                Check(ref passed, ref failed, atOnePct < 10,
                    $"Equal T/D scores → few at 1% FDR ({atOnePct})");
            }

            // Edge 3: Peptide-level q-values shared across charge states
            {
                var results = new List<DiaSearchResult>();
                for (int z = 2; z <= 4; z++)
                {
                    var r = new DiaSearchResult("HAPPYPEPTIDE", z, 500.0 / z, 0, false, 6, 30.0, 25f, 35f);
                    r.SpectralAngleScore = 0.95f - (z - 2) * 0.05f;
                    r.FragmentsDetected = 6;
                    results.Add(r);
                }
                for (int i = 0; i < 5; i++)
                {
                    var d = new DiaSearchResult($"DECOY{i}", 2, 600.0 + i, 0, true, 6, 30.0, 25f, 35f);
                    d.SpectralAngleScore = 0.3f + i * 0.05f;
                    d.FragmentsDetected = 4;
                    results.Add(d);
                }
                RunFdr(results);

                var happy = results.Where(r => r.Sequence == "HAPPYPEPTIDE").ToList();
                double pepQv = happy[0].FdrInfo.PeptideQValue!.Value;
                bool shared = happy.All(r => Math.Abs(r.FdrInfo.PeptideQValue!.Value - pepQv) < 1e-10);
                Check(ref passed, ref failed, shared,
                    $"Peptide q-value shared across charge states ({pepQv:F6})");
            }

            // Edge 4: Single result
            {
                var results = new List<DiaSearchResult>
                {
                    new DiaSearchResult("SOLO", 2, 500.0, 0, false, 6, 30.0, 25f, 35f)
                    {
                        SpectralAngleScore = 0.99f,
                        FragmentsDetected = 6
                    }
                };
                RunFdr(results);

                Check(ref passed, ref failed,
                    results[0].FdrInfo != null && results[0].FdrInfo.QValue == 0.0,
                    "Single target → q-value = 0");
            }

            // Edge 5: All decoys
            {
                var results = new List<DiaSearchResult>();
                for (int i = 0; i < 10; i++)
                {
                    var d = new DiaSearchResult($"DECOY{i}", 2, 500.0 + i, 0, true, 6, 30.0, 25f, 35f);
                    d.SpectralAngleScore = 0.5f + i * 0.03f;
                    d.FragmentsDetected = 4;
                    results.Add(d);
                }
                RunFdr(results);

                bool allHigh = results.All(r => r.FdrInfo.QValue >= 1.0);
                Check(ref passed, ref failed, allHigh,
                    "All decoys → all q-values = 1.0");
            }

            Console.WriteLine($"  Edge cases: {passed} passed, {failed} failed");
            Console.WriteLine();
        }

        // ──────────────────────────────────────────────────────────────────────
        //  FDR helper (standalone, no MetaMorpheus dependency)
        // ──────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Runs target-decoy FDR on DiaSearchResult list using spectral angle score.
        /// This is a standalone implementation that mirrors FdrAnalysisEngineDia logic
        /// so that this benchmark can run in the mzLib Development project without
        /// referencing MetaMorpheus EngineLayer.
        /// </summary>
        private static FdrSummary RunFdr(List<DiaSearchResult> results)
        {
            if (results.Count == 0)
                return new FdrSummary(0, 0, 0, 0);

            // Sort by spectral angle descending, tiebreak on FragmentsDetected
            results.Sort((a, b) =>
            {
                float sa = float.IsNaN(a.SpectralAngleScore) ? -1f : a.SpectralAngleScore;
                float sb = float.IsNaN(b.SpectralAngleScore) ? -1f : b.SpectralAngleScore;
                int cmp = sb.CompareTo(sa);
                if (cmp != 0) return cmp;
                return b.FragmentsDetected.CompareTo(a.FragmentsDetected);
            });

            // Forward pass: cumulative T/D counts and raw FDR
            int cumTarget = 0;
            int cumDecoy = 0;
            for (int i = 0; i < results.Count; i++)
            {
                if (results[i].IsDecoy) cumDecoy++;
                else cumTarget++;

                results[i].FdrInfo = new DiaFdrInfo
                {
                    CumulativeTarget = cumTarget,
                    CumulativeDecoy = cumDecoy,
                    QValue = (double)cumDecoy / Math.Max(cumTarget, 1)
                };
            }

            // Backward pass: monotonize q-values
            double runningMin = results[results.Count - 1].FdrInfo.QValue;
            for (int i = results.Count - 1; i >= 0; i--)
            {
                if (results[i].FdrInfo.QValue < runningMin)
                    runningMin = results[i].FdrInfo.QValue;
                else
                    results[i].FdrInfo.QValue = runningMin;
            }

            // Cap at 1.0
            foreach (var r in results)
                if (r.FdrInfo.QValue > 1.0) r.FdrInfo.QValue = 1.0;

            // Peptide-level: best score per sequence, recalculate, propagate
            var bestPerSequence = results
                .GroupBy(r => r.Sequence)
                .Select(g => g.OrderByDescending(r => float.IsNaN(r.SpectralAngleScore) ? -1f : r.SpectralAngleScore)
                              .First())
                .OrderByDescending(r => float.IsNaN(r.SpectralAngleScore) ? -1f : r.SpectralAngleScore)
                .ToList();

            int pepCumT = 0, pepCumD = 0;
            foreach (var r in bestPerSequence)
            {
                if (r.IsDecoy) pepCumD++; else pepCumT++;
                r.FdrInfo.PeptideQValue = (double)pepCumD / Math.Max(pepCumT, 1);
            }

            // Monotonize peptide q-values
            double pepMin = bestPerSequence[bestPerSequence.Count - 1].FdrInfo.PeptideQValue!.Value;
            for (int i = bestPerSequence.Count - 1; i >= 0; i--)
            {
                double v = bestPerSequence[i].FdrInfo.PeptideQValue!.Value;
                if (v < pepMin) pepMin = v;
                else bestPerSequence[i].FdrInfo.PeptideQValue = pepMin;
            }

            // Propagate peptide q-values
            var pepQvMap = new Dictionary<string, double>();
            foreach (var r in bestPerSequence)
            {
                double qv = Math.Min(r.FdrInfo.PeptideQValue!.Value, 1.0);
                pepQvMap[r.Sequence] = qv;
            }
            foreach (var r in results)
            {
                if (pepQvMap.TryGetValue(r.Sequence, out double pqv))
                    r.FdrInfo.PeptideQValue = pqv;
            }

            int targets1Pct = results.Count(r => !r.IsDecoy && r.FdrInfo.QValue <= 0.01);
            int peptides1Pct = pepQvMap.Values.Count(v => v <= 0.01);
            return new FdrSummary(cumTarget, cumDecoy, targets1Pct, peptides1Pct);
        }

        private readonly struct FdrSummary
        {
            public readonly int TotalTargets;
            public readonly int TotalDecoys;
            public readonly int TargetsAt1Pct;
            public readonly int PeptidesAt1Pct;

            public FdrSummary(int totalTargets, int totalDecoys, int targetsAt1Pct, int peptidesAt1Pct)
            {
                TotalTargets = totalTargets;
                TotalDecoys = totalDecoys;
                TargetsAt1Pct = targetsAt1Pct;
                PeptidesAt1Pct = peptidesAt1Pct;
            }
        }

        // ──────────────────────────────────────────────────────────────────────
        //  TSV writer (standalone, no MetaMorpheus dependency)
        // ──────────────────────────────────────────────────────────────────────

        private static string WriteTsv(List<DiaSearchResult> results, string path)
        {
            using var writer = new StreamWriter(path);

            writer.WriteLine(string.Join("\t",
                "Sequence", "Charge", "Precursor m/z", "Window ID", "Is Decoy",
                "Dot Product Score", "Spectral Angle Score",
                "Fragments Detected", "Fragments Queried", "Fragment Detection Rate",
                "Library RT", "RT Window Start", "RT Window End",
                "QValue", "Peptide QValue"));

            foreach (var r in results)
            {
                double detRate = r.FragmentsQueried > 0
                    ? (double)r.FragmentsDetected / r.FragmentsQueried
                    : 0;

                writer.WriteLine(string.Join("\t",
                    r.Sequence,
                    r.ChargeState.ToString(CultureInfo.InvariantCulture),
                    r.PrecursorMz.ToString("F4", CultureInfo.InvariantCulture),
                    r.WindowId.ToString(CultureInfo.InvariantCulture),
                    r.IsDecoy ? "True" : "False",
                    r.DotProductScore.ToString("F6", CultureInfo.InvariantCulture),
                    r.SpectralAngleScore.ToString("F6", CultureInfo.InvariantCulture),
                    r.FragmentsDetected.ToString(CultureInfo.InvariantCulture),
                    r.FragmentsQueried.ToString(CultureInfo.InvariantCulture),
                    detRate.ToString("F4", CultureInfo.InvariantCulture),
                    r.LibraryRetentionTime?.ToString("F2", CultureInfo.InvariantCulture) ?? "",
                    r.RtWindowStart.ToString("F2", CultureInfo.InvariantCulture),
                    r.RtWindowEnd.ToString("F2", CultureInfo.InvariantCulture),
                    r.FdrInfo?.QValue.ToString("F6", CultureInfo.InvariantCulture) ?? "",
                    r.FdrInfo?.PeptideQValue?.ToString("F6", CultureInfo.InvariantCulture) ?? ""));
            }

            return path;
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Synthetic data generation
        // ──────────────────────────────────────────────────────────────────────

        /// <summary>
        /// Generates synthetic DIA data with planted target signals.
        ///
        /// For each target:
        ///   - Picks a window, generates fragment m/z values
        ///   - Plants those fragments into scans near the target's RT (Gaussian falloff)
        ///
        /// For each decoy:
        ///   - Uses same structure but shifts fragment m/z by +50 Da
        ///   - These miss planted signals → decoys score poorly
        /// </summary>
        private static (MsDataScan[] Scans, List<LibraryPrecursorInput> Targets, List<LibraryPrecursorInput> Decoys)
            GenerateSyntheticDiaData(
                int windowCount, int scansPerWindow, int peaksPerScan,
                int targetCount, int decoyCount, int fragmentsPerPrecursor,
                double windowStart, double windowWidth, Random rng)
        {
            double rtMax = 60.0;
            int totalScans = windowCount * scansPerWindow;

            // Build scan framework: each window gets evenly spaced scans, interleaved
            var scanList = new List<(int windowIdx, double rt, double center, double low, double high)>(totalScans);
            for (int w = 0; w < windowCount; w++)
            {
                double wLow = windowStart + w * windowWidth;
                double wHigh = wLow + windowWidth;
                double wCenter = (wLow + wHigh) / 2.0;

                for (int s = 0; s < scansPerWindow; s++)
                {
                    double rt = (s * windowCount + w) * rtMax / totalScans;
                    scanList.Add((w, rt, wCenter, wLow, wHigh));
                }
            }
            scanList.Sort((a, b) => a.rt.CompareTo(b.rt));

            // Generate targets + plant signals
            var targets = new List<LibraryPrecursorInput>(targetCount);
            var decoys = new List<LibraryPrecursorInput>(decoyCount);
            var plantedSignals = new List<(int scanIdx, double mz, double intensity)>();

            for (int t = 0; t < targetCount; t++)
            {
                int windowIdx = t % windowCount;
                double wLow = windowStart + windowIdx * windowWidth;
                double precursorMz = wLow + windowWidth * 0.5 + (rng.NextDouble() - 0.5) * windowWidth * 0.4;
                double rt = 5.0 + rng.NextDouble() * 50.0;

                var fragMzs = new float[fragmentsPerPrecursor];
                var fragIntensities = new float[fragmentsPerPrecursor];
                for (int f = 0; f < fragmentsPerPrecursor; f++)
                {
                    fragMzs[f] = (float)(wLow + 2.0 + f * (windowWidth - 4.0) / fragmentsPerPrecursor
                        + rng.NextDouble() * 0.5);
                    fragIntensities[f] = (float)(1000.0 * (1.0 + rng.NextDouble()));
                }

                targets.Add(new LibraryPrecursorInput(
                    $"TARGET_PEP_{t}", precursorMz, 2, rt, false, fragMzs, fragIntensities));

                // Plant in nearby scans
                for (int si = 0; si < scanList.Count; si++)
                {
                    if (scanList[si].windowIdx != windowIdx) continue;
                    double rtDiff = Math.Abs(scanList[si].rt - rt);
                    if (rtDiff > 3.0) continue;

                    double rtFactor = Math.Exp(-rtDiff * rtDiff / 2.0);
                    for (int f = 0; f < fragmentsPerPrecursor; f++)
                    {
                        double intensity = fragIntensities[f] * rtFactor * (0.8 + 0.4 * rng.NextDouble());
                        if (intensity > 10)
                            plantedSignals.Add((si, fragMzs[f], intensity));
                    }
                }

                // Paired decoy
                if (t < decoyCount)
                {
                    var decoyFragMzs = new float[fragmentsPerPrecursor];
                    var decoyFragIntensities = new float[fragmentsPerPrecursor];
                    for (int f = 0; f < fragmentsPerPrecursor; f++)
                    {
                        decoyFragMzs[f] = fragMzs[f] + 50.0f;
                        decoyFragIntensities[f] = fragIntensities[f];
                    }
                    decoys.Add(new LibraryPrecursorInput(
                        $"DECOY_PEP_{t}", precursorMz + 0.5, 2,
                        rt + rng.NextDouble() * 2.0 - 1.0,
                        true, decoyFragMzs, decoyFragIntensities));
                }
            }

            // Group planted signals by scan for efficient insertion
            var signalsByScan = plantedSignals
                .GroupBy(s => s.scanIdx)
                .ToDictionary(g => g.Key, g => g.ToList());

            // Build MsDataScan[]
            var scans = new MsDataScan[scanList.Count];
            for (int i = 0; i < scanList.Count; i++)
            {
                var (_, rt, center, low, high) = scanList[i];

                var mzList = new List<double>(peaksPerScan + 20);
                var intList = new List<double>(peaksPerScan + 20);

                // Background noise
                for (int p = 0; p < peaksPerScan; p++)
                {
                    mzList.Add(low + rng.NextDouble() * windowWidth);
                    intList.Add(rng.NextDouble() * 100.0);
                }

                // Planted signals
                if (signalsByScan.TryGetValue(i, out var signals))
                {
                    foreach (var (_, mz, intensity) in signals)
                    {
                        mzList.Add(mz);
                        intList.Add(intensity);
                    }
                }

                // Sort by m/z
                var pairs = mzList.Zip(intList, (m, inten) => (m, inten)).OrderBy(x => x.m).ToList();
                double[] mzArr = pairs.Select(x => x.m).ToArray();
                double[] intArr = pairs.Select(x => x.inten).ToArray();

                scans[i] = new MsDataScan(
                    massSpectrum: new MzSpectrum(mzArr, intArr, false),
                    oneBasedScanNumber: i + 1,
                    msnOrder: 2,
                    isCentroid: true,
                    polarity: Polarity.Positive,
                    retentionTime: rt,
                    scanWindowRange: new MzRange(low, high),
                    scanFilter: "FTMS",
                    mzAnalyzer: MZAnalyzerType.Orbitrap,
                    totalIonCurrent: intArr.Sum(),
                    injectionTime: 20.0,
                    noiseData: null,
                    nativeId: $"scan={i + 1}",
                    isolationMZ: center,
                    isolationWidth: windowWidth,
                    dissociationType: DissociationType.HCD);
            }

            return (scans, targets, decoys);
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Helpers
        // ──────────────────────────────────────────────────────────────────────

        private static void Check(ref int passed, ref int failed, bool condition, string description)
        {
            if (condition)
            {
                Console.WriteLine($"  ✓ {description}");
                passed++;
            }
            else
            {
                Console.WriteLine($"  ✗ FAIL: {description}");
                failed++;
            }
        }
    }
}