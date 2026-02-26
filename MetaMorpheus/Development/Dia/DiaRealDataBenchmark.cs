// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using EngineLayer;
using EngineLayer.DiaSearch;
using EngineLayer.FdrAnalysisDia;
using MassSpectrometry;
using MassSpectrometry.Dia;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;

namespace Development.Dia
{
    /// <summary>
    /// End-to-end benchmark using real DIA data files:
    ///   - E24484_Mag_1_4_50_3.mzML  (~20 MB, 767 scans, 152 DIA windows, RT 37–38 min)
    ///   - target.msp                 (2994 target library entries, avg ~21 fragments)
    ///   - decoy.msp                  (3020 decoy library entries, avg ~21 fragments)
    ///
    /// Exercises the full MetaMorpheus DIA pipeline:
    ///   mzML load → .msp parse → DiaEngine (index → query → extract → score) → FDR → TSV
    ///
    /// Data files should be placed in: DiaDevBenchmarks/TestData/
    ///   (or set MMDIADEV_TESTDATA environment variable to their folder)
    ///
    /// Run: dotnet run -c Release --project MetaMorpheus/DiaDevBenchmarks
    /// Use Ctrl+F5 (no debugger) for accurate timing.
    /// </summary>
    public static class DiaRealDataBenchmark
    {
        private const string MzmlFileName = "E24484_Mag_1_4_50_3.mzML";
        private const string TargetMspFileName = "target.msp";
        private const string DecoyMspFileName = "decoy.msp";

        public static void RunAll()
        {
            Console.WriteLine("========================================================");
            Console.WriteLine("  DIA Real-Data End-to-End Benchmark");
            Console.WriteLine("========================================================");
            Console.WriteLine();
            Console.WriteLine($"  Processors:    {Environment.ProcessorCount}");
            Console.WriteLine($"  OS:            {RuntimeInformation.OSDescription}");
            Console.WriteLine($"  Runtime:       {RuntimeInformation.FrameworkDescription}");
            Console.WriteLine($"  GPU backend:   {FragmentExtractorFactory.DescribeBackend()}");
            Console.WriteLine();

            string dataDir = ResolveDataDirectory();
            if (dataDir == null) return;

            string mzmlPath = Path.Combine(dataDir, MzmlFileName);
            string targetMspPath = Path.Combine(dataDir, TargetMspFileName);
            string decoyMspPath = Path.Combine(dataDir, DecoyMspFileName);

            foreach (var (path, label) in new[]
            {
                (mzmlPath, "mzML"),
                (targetMspPath, "target .msp"),
                (decoyMspPath, "decoy .msp")
            })
            {
                if (!File.Exists(path))
                {
                    Console.Error.WriteLine($"  [SKIP] {label} not found: {path}");
                    Console.Error.WriteLine("  Place data files in DiaDevBenchmarks/TestData/");
                    Console.Error.WriteLine("  or set MMDIADEV_TESTDATA environment variable.");
                    return;
                }
            }

            Console.WriteLine($"  mzML file:     {Path.GetFileName(mzmlPath)}");
            Console.WriteLine($"  Size:          {new FileInfo(mzmlPath).Length / (1024.0 * 1024.0):F1} MB");
            Console.WriteLine($"  Target lib:    {Path.GetFileName(targetMspPath)}");
            Console.WriteLine($"  Decoy lib:     {Path.GetFileName(decoyMspPath)}");
            Console.WriteLine();

            // ── Load raw data ─────────────────────────────────────────────
            var loadSw = Stopwatch.StartNew();
            MsDataFile dataFile = MsDataFileReader.GetDataFile(mzmlPath);
            dataFile.LoadAllStaticData();
            MsDataScan[] scans = dataFile.GetAllScansList().ToArray();
            loadSw.Stop();

            int ms1Count = scans.Count(s => s.MsnOrder == 1);
            int ms2Count = scans.Count(s => s.MsnOrder == 2);
            Console.WriteLine($"  mzML load:     {loadSw.ElapsedMilliseconds,6} ms  ({scans.Length} scans: {ms1Count} MS1, {ms2Count} MS2)");

            // ── Parse .msp libraries ──────────────────────────────────────
            var parseSw = Stopwatch.StartNew();
            var targets = MspParser.Parse(targetMspPath, isDecoy: false);
            var decoys = MspParser.Parse(decoyMspPath, isDecoy: true);
            parseSw.Stop();

            int totalFragments = targets.Sum(t => t.FragmentMzs.Length) + decoys.Sum(d => d.FragmentMzs.Length);
            Console.WriteLine($"  .msp parse:    {parseSw.ElapsedMilliseconds,6} ms  " +
                $"({targets.Count} T + {decoys.Count} D = {targets.Count + decoys.Count} entries, " +
                $"{totalFragments:N0} fragments)");
            Console.WriteLine();

            var allLibraryEntries = new List<MspEntry>(targets.Count + decoys.Count);
            allLibraryEntries.AddRange(targets);
            allLibraryEntries.AddRange(decoys);

            // ── Run benchmarks ────────────────────────────────────────────
            RunDiaEnginePipeline(scans, allLibraryEntries);
            Console.WriteLine(new string('─', 60));
            RunMzLibLayerPipeline(scans, allLibraryEntries);
            Console.WriteLine(new string('─', 60));
            RunThreadScaling(scans, allLibraryEntries);
            Console.WriteLine(new string('─', 60));
            RunValidation(scans, allLibraryEntries);

            Console.WriteLine();
            Console.WriteLine("Done.");
        }

        // ──────────────────────────────────────────────────────────────────
        //  Benchmark 1: Full DiaEngine pipeline (MetaMorpheus facade)
        // ──────────────────────────────────────────────────────────────────

        private static void RunDiaEnginePipeline(MsDataScan[] scans, List<MspEntry> library)
        {
            Console.WriteLine("── DiaEngine Full Pipeline ─────────────────────────────────");

            var librarySpectra = ConvertToLibrarySpectra(library);
            var mzLibParams = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 5.0f,
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0.0f,
                MaxThreads = -1,
                PreferGpu = false
            };

            var commonParams = new CommonParameters();
            var nestedIds = new List<string> { "DiaDevBenchmark" };

            var totalSw = Stopwatch.StartNew();

            using var diaEngine = new DiaEngine(
                scans, librarySpectra, mzLibParams, commonParams, null, nestedIds);

            var engineResults = (DiaEngineResults)diaEngine.Run();
            totalSw.Stop();

            var results = engineResults.DiaResults;
            int targetCount = results.Count(r => !r.IsDecoy);
            int decoyCount = results.Count(r => r.IsDecoy);

            Console.WriteLine($"  DiaEngine:     {totalSw.ElapsedMilliseconds,6} ms  " +
                $"({targetCount} T, {decoyCount} D)");

            // FDR
            var fdrSw = Stopwatch.StartNew();
            var fdrEngine = new FdrAnalysisEngineDia(
                results, DiaFdrScoreType.SpectralAngle,
                commonParams, null, nestedIds);
            var fdrResults = (FdrAnalysisResultsDia)fdrEngine.Run();
            fdrSw.Stop();

            Console.WriteLine($"  FDR engine:    {fdrSw.ElapsedMilliseconds,6} ms  " +
                $"({fdrResults.ResultsWithin1PercentFdr} @ 1%, " +
                $"{fdrResults.PeptidesWithin1PercentFdr} peptides @ 1%)");

            // TSV
            string outputDir = Path.Combine(Path.GetTempPath(), "DiaDevBench_" + Guid.NewGuid().ToString("N")[..8]);
            Directory.CreateDirectory(outputDir);
            var tsvSw = Stopwatch.StartNew();
            WriteBenchmarkTsv(results, Path.Combine(outputDir, "AllDiaResults.tsv"));
            tsvSw.Stop();
            int tsvRows = File.ReadAllLines(Path.Combine(outputDir, "AllDiaResults.tsv")).Length - 1;
            Console.WriteLine($"  TSV write:     {tsvSw.ElapsedMilliseconds,6} ms  ({tsvRows} rows)");

            long totalMs = totalSw.ElapsedMilliseconds + fdrSw.ElapsedMilliseconds + tsvSw.ElapsedMilliseconds;
            Console.WriteLine($"  TOTAL:         {totalMs,6} ms");

            // Score stats
            PrintScoreStats(results);
            Console.WriteLine();

            try { Directory.Delete(outputDir, true); } catch { }
        }

        // ──────────────────────────────────────────────────────────────────
        //  Benchmark 2: mzLib layer-by-layer pipeline (detailed timing)
        // ──────────────────────────────────────────────────────────────────

        private static void RunMzLibLayerPipeline(MsDataScan[] scans, List<MspEntry> library)
        {
            Console.WriteLine("── mzLib Layer-by-Layer Pipeline ───────────────────────────");

            var precursors = ConvertToLibraryPrecursorInputs(library);
            var parameters = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 5.0f,
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0.0f,
                MaxThreads = -1,
                PreferGpu = false
            };

            var totalSw = Stopwatch.StartNew();

            // 1. Index build
            var sw = Stopwatch.StartNew();
            using var index = DiaScanIndexBuilder.Build(scans);
            sw.Stop();
            long totalPeaks = 0;
            for (int i = 0; i < index.ScanCount; i++)
                totalPeaks += index.GetScanMzSpan(i).Length;
            Console.WriteLine($"  Index build:     {sw.ElapsedMilliseconds,6} ms  " +
                $"({index.WindowCount} win, {index.ScanCount} scans, {totalPeaks:N0} peaks)");

            // 2. Query generation
            sw.Restart();
            var genResult = DiaLibraryQueryGenerator.Generate(precursors, index, parameters);
            sw.Stop();
            Console.WriteLine($"  Query gen:       {sw.ElapsedMilliseconds,6} ms  " +
                $"({genResult.Queries.Length:N0} queries, {genResult.PrecursorGroups.Length:N0} groups, " +
                $"{genResult.SkippedNoWindow} skipped)");

            // 3. Extraction
            sw.Restart();
            using var orchestrator = new DiaExtractionOrchestrator(index);
            var extraction = orchestrator.ExtractAll(genResult.Queries);
            sw.Stop();
            double qps = genResult.Queries.Length / Math.Max(0.001, sw.Elapsed.TotalSeconds);
            Console.WriteLine($"  Extraction:      {sw.ElapsedMilliseconds,6} ms  " +
                $"({extraction.TotalDataPoints:N0} data points, {qps:N0} queries/sec)");

            // 4. Score + assemble
            sw.Restart();
            var dotScorer = new NormalizedDotProductScorer();
            var saScorer = new SpectralAngleScorer();
            var results = DiaLibraryQueryGenerator.AssembleResults(
                precursors, genResult, extraction.Results, parameters,
                dotProductScorer: dotScorer, spectralAngleScorer: saScorer);
            sw.Stop();

            int targetResults = results.Count(r => !r.IsDecoy);
            int decoyResults = results.Count(r => r.IsDecoy);
            Console.WriteLine($"  Score+assemble:  {sw.ElapsedMilliseconds,6} ms  " +
                $"({results.Count} results: {targetResults} T, {decoyResults} D)");

            totalSw.Stop();
            Console.WriteLine($"  TOTAL:           {totalSw.ElapsedMilliseconds,6} ms");
            Console.WriteLine();
        }

        // ──────────────────────────────────────────────────────────────────
        //  Benchmark 3: Thread scaling
        // ──────────────────────────────────────────────────────────────────

        private static void RunThreadScaling(MsDataScan[] scans, List<MspEntry> library)
        {
            Console.WriteLine("── Thread Scaling ─────────────────────────────────────────");

            var precursors = ConvertToLibraryPrecursorInputs(library);
            var parameters = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 5.0f,
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0.0f,
                PreferGpu = false
            };

            using var index = DiaScanIndexBuilder.Build(scans);
            var genResult = DiaLibraryQueryGenerator.Generate(precursors, index, parameters);

            int[] threadCounts = { 1, 2, 4, 8, 16, 32 };
            threadCounts = threadCounts.Where(t => t <= Environment.ProcessorCount * 2).ToArray();

            double baselineMs = 0;
            foreach (int threads in threadCounts)
            {
                // Warmup
                using (var warmup = new DiaExtractionOrchestrator(index))
                    warmup.ExtractAll(genResult.Queries, maxDegreeOfParallelism: threads);

                // Timed (avg of 3)
                const int runs = 3;
                double totalMs = 0;
                for (int r = 0; r < runs; r++)
                {
                    var sw = Stopwatch.StartNew();
                    using var orch = new DiaExtractionOrchestrator(index);
                    orch.ExtractAll(genResult.Queries, maxDegreeOfParallelism: threads);
                    sw.Stop();
                    totalMs += sw.Elapsed.TotalMilliseconds;
                }

                double avgMs = totalMs / runs;
                if (threads == 1) baselineMs = avgMs;
                double speedup = baselineMs / Math.Max(avgMs, 0.01);

                Console.WriteLine($"  {threads,2} threads:  {avgMs,8:F1} ms  |  {speedup,5:F2}x speedup");
            }
            Console.WriteLine();
        }

        // ──────────────────────────────────────────────────────────────────
        //  Benchmark 4: Validation checks
        // ──────────────────────────────────────────────────────────────────

        private static void RunValidation(MsDataScan[] scans, List<MspEntry> library)
        {
            Console.WriteLine("── Validation ─────────────────────────────────────────────");
            int passed = 0, failed = 0;

            var librarySpectra = ConvertToLibrarySpectra(library);
            var mzLibParams = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 5.0f,
                MinFragmentsRequired = 3,
                MinScoreThreshold = 0.0f,
                MaxThreads = -1,
                PreferGpu = false
            };

            // Run DiaEngine
            var commonParams = new CommonParameters();
            var nestedIds = new List<string> { "Validation" };

            using var diaEngine = new DiaEngine(
                scans, librarySpectra, mzLibParams, commonParams, null, nestedIds);
            var engineResults = (DiaEngineResults)diaEngine.Run();
            var results = engineResults.DiaResults;

            // Run FDR
            var fdrEngine = new FdrAnalysisEngineDia(
                results, DiaFdrScoreType.SpectralAngle,
                commonParams, null, nestedIds);
            var fdrResults = (FdrAnalysisResultsDia)fdrEngine.Run();

            var targetResults = results.Where(r => !r.IsDecoy).ToList();
            var decoyResults = results.Where(r => r.IsDecoy).ToList();

            // 1. Both T and D produced
            Check(ref passed, ref failed,
                targetResults.Count > 0 && decoyResults.Count > 0,
                $"Both T and D have results ({targetResults.Count} T, {decoyResults.Count} D)");

            // 2. Targets score higher than decoys on average
            if (targetResults.Count > 0 && decoyResults.Count > 0)
            {
                double avgT = targetResults.Where(r => !float.IsNaN(r.SpectralAngleScore))
                    .Average(r => (double)r.SpectralAngleScore);
                double avgD = decoyResults.Where(r => !float.IsNaN(r.SpectralAngleScore))
                    .Average(r => (double)r.SpectralAngleScore);
                Check(ref passed, ref failed, avgT > avgD,
                    $"Target avg SA ({avgT:F4}) > Decoy avg SA ({avgD:F4})");
            }

            // 3. Q-values monotonically non-decreasing
            var sorted = results
                .OrderByDescending(r => float.IsNaN(r.SpectralAngleScore) ? -1f : r.SpectralAngleScore)
                .ToList();
            bool monotonic = true;
            for (int i = 1; i < sorted.Count; i++)
            {
                if (sorted[i].FdrInfo != null && sorted[i - 1].FdrInfo != null &&
                    sorted[i].FdrInfo.QValue < sorted[i - 1].FdrInfo.QValue - 1e-10)
                {
                    monotonic = false;
                    break;
                }
            }
            Check(ref passed, ref failed, monotonic, "Q-values monotonically non-decreasing");

            // 4. Some targets pass 1% FDR
            Check(ref passed, ref failed,
                fdrResults.ResultsWithin1PercentFdr > 0,
                $"Targets at 1% FDR: {fdrResults.ResultsWithin1PercentFdr}");

            // 5. Peptides at 1% FDR
            Check(ref passed, ref failed,
                fdrResults.PeptidesWithin1PercentFdr > 0,
                $"Peptides at 1% FDR: {fdrResults.PeptidesWithin1PercentFdr}");

            // 6. All results have FdrInfo
            Check(ref passed, ref failed,
                results.All(r => r.FdrInfo != null),
                "All results have FdrInfo");

            // 7. All results have peptide q-values
            Check(ref passed, ref failed,
                results.All(r => r.FdrInfo?.PeptideQValue != null),
                "All results have peptide q-values");

            // 8. Q-values in valid range
            Check(ref passed, ref failed,
                results.Where(r => r.FdrInfo != null)
                    .All(r => r.FdrInfo.QValue >= 0 && r.FdrInfo.QValue <= 1.0),
                "All q-values in [0, 1]");

            // 9. DiaEngine target/decoy counts match
            Check(ref passed, ref failed,
                engineResults.TargetCount == targetResults.Count &&
                engineResults.DecoyCount == decoyResults.Count,
                $"DiaEngineResults counts match (T={engineResults.TargetCount}, D={engineResults.DecoyCount})");

            // 10. TSV round-trip
            string tmpDir = Path.Combine(Path.GetTempPath(), "DiaVal_" + Guid.NewGuid().ToString("N")[..8]);
            Directory.CreateDirectory(tmpDir);
            string tsvPath = WriteBenchmarkTsv(results, Path.Combine(tmpDir, "test.tsv"));
            var lines = File.ReadAllLines(tsvPath);
            string[] hdCols = lines[0].Split('\t');
            int qvIdx = Array.IndexOf(hdCols, "QValue");
            bool parseable = lines.Length > 1 && qvIdx >= 0;
            if (parseable)
            {
                string[] row = lines[1].Split('\t');
                parseable = row.Length == hdCols.Length &&
                    double.TryParse(row[qvIdx], NumberStyles.Float, CultureInfo.InvariantCulture, out _);
            }
            Check(ref passed, ref failed, parseable, "TSV round-trip parseable");
            try { Directory.Delete(tmpDir, true); } catch { }

            Console.WriteLine();
            Console.WriteLine($"  Validation: {passed} passed, {failed} failed");
            Console.WriteLine();
        }

        // ──────────────────────────────────────────────────────────────────
        //  Helpers
        // ──────────────────────────────────────────────────────────────────

        /// <summary>
        /// Converts parsed .msp entries to LibrarySpectrum objects (what DiaEngine expects).
        /// Creates MatchedFragmentIon objects for each fragment.
        /// </summary>
        private static List<LibrarySpectrum> ConvertToLibrarySpectra(List<MspEntry> entries)
        {
            var result = new List<LibrarySpectrum>(entries.Count);
            foreach (var e in entries)
            {
                var matchedIons = new List<MatchedFragmentIon>(e.FragmentMzs.Length);
                for (int i = 0; i < e.FragmentMzs.Length; i++)
                {
                    // Create a minimal MatchedFragmentIon with the m/z and intensity
                    // The Product fields don't matter for DIA extraction — only Mz and Intensity are used
                    var product = new Product(ProductType.y, FragmentationTerminus.C, 0, 0, 0, 0);
                    matchedIons.Add(new MatchedFragmentIon(product, e.FragmentMzs[i], e.FragmentIntensities[i], 1));
                }

                var spectrum = new LibrarySpectrum(
                    e.Sequence, e.PrecursorMz, e.Charge, matchedIons, e.RetentionTime, e.IsDecoy);
                result.Add(spectrum);
            }
            return result;
        }

        /// <summary>
        /// Converts parsed .msp entries to the LibraryPrecursorInput format expected by mzLib DIA engine directly.
        /// </summary>
        private static List<LibraryPrecursorInput> ConvertToLibraryPrecursorInputs(List<MspEntry> entries)
        {
            return entries.Select(e => new LibraryPrecursorInput(
                e.Sequence, e.PrecursorMz, e.Charge, e.RetentionTime,
                e.IsDecoy, e.FragmentMzs, e.FragmentIntensities)).ToList();
        }

        private static void PrintScoreStats(List<DiaSearchResult> results)
        {
            var targetScores = results.Where(r => !r.IsDecoy && !float.IsNaN(r.SpectralAngleScore))
                .Select(r => r.SpectralAngleScore).ToList();
            var decoyScores = results.Where(r => r.IsDecoy && !float.IsNaN(r.SpectralAngleScore))
                .Select(r => r.SpectralAngleScore).ToList();

            if (targetScores.Count > 0)
            {
                targetScores.Sort();
                Console.WriteLine($"  Target SA:     min={targetScores[0]:F4}  " +
                    $"med={targetScores[targetScores.Count / 2]:F4}  " +
                    $"max={targetScores[^1]:F4}  (n={targetScores.Count})");
            }
            if (decoyScores.Count > 0)
            {
                decoyScores.Sort();
                Console.WriteLine($"  Decoy SA:      min={decoyScores[0]:F4}  " +
                    $"med={decoyScores[decoyScores.Count / 2]:F4}  " +
                    $"max={decoyScores[^1]:F4}  (n={decoyScores.Count})");
            }
        }

        private static string WriteBenchmarkTsv(List<DiaSearchResult> results, string path)
        {
            using var writer = new StreamWriter(path);
            writer.WriteLine(string.Join("\t",
                "Sequence", "Charge", "Precursor m/z", "Window ID", "Is Decoy",
                "Dot Product Score", "Spectral Angle Score",
                "Fragments Detected", "Fragments Queried", "Fragment Detection Rate",
                "Library RT", "RT Window Start", "RT Window End",
                "QValue", "Peptide QValue"));

            foreach (var r in results.OrderBy(r => r.FdrInfo?.QValue ?? 2.0))
            {
                writer.WriteLine(string.Join("\t",
                    r.Sequence,
                    r.ChargeState.ToString(CultureInfo.InvariantCulture),
                    r.PrecursorMz.ToString("F4", CultureInfo.InvariantCulture),
                    r.WindowId.ToString(CultureInfo.InvariantCulture),
                    r.IsDecoy ? "True" : "False",
                    float.IsNaN(r.DotProductScore) ? "" : r.DotProductScore.ToString("F6", CultureInfo.InvariantCulture),
                    float.IsNaN(r.SpectralAngleScore) ? "" : r.SpectralAngleScore.ToString("F6", CultureInfo.InvariantCulture),
                    r.FragmentsDetected.ToString(CultureInfo.InvariantCulture),
                    r.FragmentsQueried.ToString(CultureInfo.InvariantCulture),
                    r.FragmentDetectionRate.ToString("F4", CultureInfo.InvariantCulture),
                    r.LibraryRetentionTime?.ToString("F2", CultureInfo.InvariantCulture) ?? "",
                    r.RtWindowStart.ToString("F2", CultureInfo.InvariantCulture),
                    r.RtWindowEnd.ToString("F2", CultureInfo.InvariantCulture),
                    r.FdrInfo?.QValue.ToString("F6", CultureInfo.InvariantCulture) ?? "",
                    r.FdrInfo?.PeptideQValue?.ToString("F6", CultureInfo.InvariantCulture) ?? ""));
            }
            return path;
        }

        private static string ResolveDataDirectory()
        {
            string envDir = Environment.GetEnvironmentVariable("MMDIADEV_TESTDATA");
            if (!string.IsNullOrEmpty(envDir) && Directory.Exists(envDir))
                return envDir;

            // Walk up from working directory looking for Test/DiaTestData
            string[] candidates =
            {
                @"Test\DiaTestData",
                @"..\Test\DiaTestData",
                @"..\..\Test\DiaTestData",
                @"..\..\..\Test\DiaTestData",
                @"..\..\..\..\Test\DiaTestData",
                @"..\..\..\..\..\Test\DiaTestData",
            };

            foreach (string candidate in candidates)
            {
                string fullPath = Path.GetFullPath(candidate);
                if (Directory.Exists(fullPath) &&
                    File.Exists(Path.Combine(fullPath, MzmlFileName)))
                    return fullPath;
            }

            Console.Error.WriteLine("  [SKIP] Test data directory not found.");
            Console.Error.WriteLine(@"  Place data files in Test\DiaTestData\:");
            Console.Error.WriteLine($"    - {MzmlFileName}");
            Console.Error.WriteLine($"    - {TargetMspFileName}");
            Console.Error.WriteLine($"    - {DecoyMspFileName}");
            Console.Error.WriteLine("  Or set MMDIADEV_TESTDATA environment variable.");
            return null;
        }

        private static void Check(ref int passed, ref int failed, bool condition, string description)
        {
            if (condition) { Console.WriteLine($"  ✓ {description}"); passed++; }
            else { Console.WriteLine($"  ✗ FAIL: {description}"); failed++; }
        }
    }
}