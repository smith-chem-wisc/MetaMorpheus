// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DiaSearch;
using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Readers;

namespace Development.Dia
{
    /// <summary>
    /// End-to-end benchmarks for the DiaEngine facade using real DIA mzML data.
    /// 
    /// Uses the small real DIA file E24484_Mag_1_4_50_3.mzML from mzLib test data.
    /// This file contains real DIA acquisition data with isolation windows, providing
    /// a realistic test of the full pipeline:
    ///   mzML → DiaScanIndex → QueryGeneration → Extraction → Scoring → Results
    /// </summary>
    public static class DiaEngineBenchmark
    {
        private static readonly CommonParameters DefaultCommonParameters = new CommonParameters();
        private static readonly List<(string, CommonParameters)> NoFileSpecificParams = new();
        private static readonly List<string> BenchmarkNestedIds = new() { "Benchmark", "DiaEngine" };

        private static string FindMzmlPath()
        {
            // Layout:  E:\GitClones\MetaMorpheus\MetaMorpheus\Development\Dia\bin\Debug\net8.0\  (CWD)
            //          E:\GitClones\mzLib\mzLib\Test\DataFiles\E24484_Mag_1_4_50_3.mzML
            // So from bin\Debug\net8.0: ..\..\..\..\..\..\..\ = E:\GitClones\
            string baseDir = AppDomain.CurrentDomain.BaseDirectory;
            string[] candidatePaths =
            {
                // From bin output: up 7 levels to GitClones root, then into mzLib\mzLib\Test\DataFiles
                Path.Combine(baseDir, "..", "..", "..", "..", "..", "..", "..", "mzLib", "mzLib", "Test", "DataFiles", "E24484_Mag_1_4_50_3.mzML"),
                // From bin output: up 7 levels, mzLib without nested mzLib (flat clone)
                Path.Combine(baseDir, "..", "..", "..", "..", "..", "..", "..", "mzLib", "Test", "DataFiles", "E24484_Mag_1_4_50_3.mzML"),
                // From project dir (dotnet run): up 4 levels to GitClones root
                Path.Combine(baseDir, "..", "..", "..", "..", "mzLib", "mzLib", "Test", "DataFiles", "E24484_Mag_1_4_50_3.mzML"),
                Path.Combine(baseDir, "..", "..", "..", "..", "mzLib", "Test", "DataFiles", "E24484_Mag_1_4_50_3.mzML"),
                // Direct environment variable override
                Path.Combine(Environment.GetEnvironmentVariable("MZLIB_TEST_DATA") ?? "", "E24484_Mag_1_4_50_3.mzML"),
            };

            foreach (var path in candidatePaths)
            {
                string fullPath = Path.GetFullPath(path);
                if (File.Exists(fullPath))
                    return fullPath;
            }
            return null;
        }

        public static void RunAll()
        {
            Console.WriteLine("========================================================");
            Console.WriteLine("  DiaEngine End-to-End Benchmark");
            Console.WriteLine("========================================================");
            Console.WriteLine();
            Console.WriteLine($"  Processors:    {Environment.ProcessorCount}");
            Console.WriteLine($"  OS:            {System.Runtime.InteropServices.RuntimeInformation.OSDescription}");
            Console.WriteLine($"  GPU backend:   {FragmentExtractorFactory.DescribeBackend()}");
            Console.WriteLine();

            string mzmlPath = FindMzmlPath();
            if (mzmlPath == null)
            {
                Console.WriteLine("  ERROR: Could not find E24484_Mag_1_4_50_3.mzML");
                Console.WriteLine("  Base directory: " + AppDomain.CurrentDomain.BaseDirectory);
                Console.WriteLine("  Tip: set MZLIB_TEST_DATA env var to the mzLib Test/DataFiles folder");
                Console.WriteLine();
                Console.WriteLine("  Falling back to synthetic-only benchmarks...");
                Console.WriteLine();
                BenchmarkSyntheticFullPipeline();
                return;
            }

            Console.WriteLine($"  mzML file:     {Path.GetFileName(mzmlPath)}");
            Console.WriteLine($"  Size:          {new FileInfo(mzmlPath).Length / 1024.0 / 1024.0:F1} MB");
            Console.WriteLine();

            var loadTimer = Stopwatch.StartNew();
            var mzml = new Mzml(mzmlPath);
            mzml.LoadAllStaticData();
            loadTimer.Stop();

            var allScans = mzml.GetMsDataScans();
            int ms1Count = allScans.Count(s => s != null && s.MsnOrder == 1);
            int ms2Count = allScans.Count(s => s != null && s.MsnOrder == 2);

            Console.WriteLine($"  Load time:     {loadTimer.ElapsedMilliseconds} ms");
            Console.WriteLine($"  Total scans:   {allScans.Length} (MS1: {ms1Count}, MS2: {ms2Count})");
            Console.WriteLine();

            BenchmarkRealDataIndexBuild(allScans);
            BenchmarkRealDataFullPipeline(allScans, precursorCount: 100);
            BenchmarkRealDataFullPipeline(allScans, precursorCount: 1_000);
            BenchmarkRealDataFullPipeline(allScans, precursorCount: 10_000);
            BenchmarkRealDataScaling(allScans);
            BenchmarkSyntheticFullPipeline();

            Console.WriteLine("Done.");
        }

        // ──────────────────────────────────────────────────────────────────────

        private static void BenchmarkRealDataIndexBuild(MsDataScan[] scans)
        {
            Console.WriteLine("── Real Data: Index Build ────────────────────────────────────");

            using (DiaScanIndexBuilder.Build(scans)) { } // warmup

            const int iterations = 5;
            var sw = Stopwatch.StartNew();
            DiaScanIndex index = null;
            for (int i = 0; i < iterations; i++)
            {
                index?.Dispose();
                index = DiaScanIndexBuilder.Build(scans);
            }
            sw.Stop();

            double avgMs = sw.Elapsed.TotalMilliseconds / iterations;

            Console.WriteLine($"  Build time:    {avgMs:F1} ms (avg of {iterations})");
            Console.WriteLine($"  Scans:         {index.ScanCount}");
            Console.WriteLine($"  Windows:       {index.WindowCount}");
            Console.WriteLine($"  Total peaks:   {index.TotalPeakCount:N0}");
            Console.WriteLine($"  Throughput:    {index.TotalPeakCount / (avgMs / 1000.0):N0} peaks/sec");
            Console.WriteLine();

            Console.WriteLine("  Window layout:");
            foreach (int wid in index.GetWindowIds().OrderBy(w => w))
            {
                var bounds = index.GetWindowBounds(wid);
                index.TryGetScanRangeForWindow(wid, out _, out int count);
                Console.WriteLine($"    Window {wid,2}: [{bounds.LowerBound:F1} – {bounds.UpperBound:F1}] m/z, {count} scans");
            }
            Console.WriteLine();

            index.Dispose();
        }

        // ──────────────────────────────────────────────────────────────────────

        private static void BenchmarkRealDataFullPipeline(MsDataScan[] scans, int precursorCount)
        {
            Console.WriteLine($"── Real Data: Full Pipeline ({precursorCount:N0} precursors) ──────────────");

            using var index = DiaScanIndexBuilder.Build(scans);
            var library = GenerateSyntheticLibrary(index, precursorCount, fragmentsPerPrecursor: 12);
            Console.WriteLine($"  Library size:  {library.Count} entries, {library.Sum(l => l.MatchedFragmentIons.Count)} total fragments");

            // Single-threaded
            RunAndPrint(scans, library, maxThreads: 1);

            // Multi-threaded
            int threads = Math.Min(Environment.ProcessorCount, 16);
            var result = RunAndPrint(scans, library, maxThreads: threads);

            Console.WriteLine($"  Matches:       {result.DiaResults.Count:N0}");
            if (result.DiaResults.Count > 0)
            {
                float avgDot = result.DiaResults
                    .Where(r => !float.IsNaN(r.DotProductScore))
                    .Select(r => r.DotProductScore)
                    .DefaultIfEmpty(0).Average();
                Console.WriteLine($"  Avg dot prod:  {avgDot:F4}");
            }
            Console.WriteLine();
        }

        // ──────────────────────────────────────────────────────────────────────

        private static void BenchmarkRealDataScaling(MsDataScan[] scans)
        {
            Console.WriteLine("── Real Data: Thread Scaling (1000 precursors) ───────────────");

            using var index = DiaScanIndexBuilder.Build(scans);
            var library = GenerateSyntheticLibrary(index, 1000, 12);

            // Warmup
            RunEngine(scans, library, maxThreads: 1);

            int[] threadCounts = { 1, 2, 4, 8, Math.Min(16, Environment.ProcessorCount) };
            double serialMs = 0;

            foreach (int threads in threadCounts)
            {
                const int iterations = 3;
                double totalMs = 0;
                DiaEngineResults lastResult = null;

                for (int i = 0; i < iterations; i++)
                {
                    var sw = Stopwatch.StartNew();
                    lastResult = RunEngine(scans, library, threads);
                    sw.Stop();
                    totalMs += sw.Elapsed.TotalMilliseconds;
                }

                double avgMs = totalMs / iterations;
                if (threads == 1) serialMs = avgMs;
                double speedup = serialMs / avgMs;

                Console.WriteLine(
                    $"  {threads,2} threads:  {avgMs,8:F1} ms  |  " +
                    $"{speedup,5:F2}× speedup  |  " +
                    $"{lastResult.DiaResults.Count} matches");
            }
            Console.WriteLine();
        }

        // ──────────────────────────────────────────────────────────────────────

        private static void BenchmarkSyntheticFullPipeline()
        {
            Console.WriteLine("── Synthetic Data: Full Pipeline ─────────────────────────────");

            int windowCount = 20, scansPerWindow = 200, peaksPerScan = 150;
            var scans = GenerateSyntheticScans(windowCount, scansPerWindow, peaksPerScan);
            Console.WriteLine($"  Scans generated: {scans.Length}");

            using var index = DiaScanIndexBuilder.Build(scans);
            var library = GenerateSyntheticLibrary(index, 5_000, 10);
            Console.WriteLine($"  Library: {library.Count} precursors");

            foreach (int pc in new[] { 100, 500, 1_000, 5_000 })
            {
                var subLib = library.GetRange(0, Math.Min(pc, library.Count));
                int threads = Math.Min(Environment.ProcessorCount, 8);

                var sw = Stopwatch.StartNew();
                var result = RunEngine(scans, subLib, threads);
                sw.Stop();

                Console.WriteLine(
                    $"  {pc,5} precursors:  {sw.ElapsedMilliseconds,5} ms  |  " +
                    $"{result.DiaResults.Count,5} matches");
            }
            Console.WriteLine();
        }

        // ──────────────────────────────────────────────────────────────────────
        //  Helpers
        // ──────────────────────────────────────────────────────────────────────

        private static DiaEngineResults RunEngine(MsDataScan[] scans, List<LibrarySpectrum> library, int maxThreads)
        {
            var diaParams = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 5f,
                MinFragmentsRequired = 3,
                MaxThreads = maxThreads
            };
            var engine = new DiaEngine(scans, library, diaParams,
                DefaultCommonParameters, NoFileSpecificParams, BenchmarkNestedIds);
            return (DiaEngineResults)engine.Run();
        }

        private static DiaEngineResults RunAndPrint(MsDataScan[] scans, List<LibrarySpectrum> library, int maxThreads)
        {
            var sw = Stopwatch.StartNew();
            var result = RunEngine(scans, library, maxThreads);
            sw.Stop();
            Console.WriteLine($"  [{maxThreads} thread{(maxThreads > 1 ? "s" : " ")}]  {sw.ElapsedMilliseconds} ms");
            return result;
        }

        private static List<LibrarySpectrum> GenerateSyntheticLibrary(
            DiaScanIndex index, int precursorCount, int fragmentsPerPrecursor)
        {
            var rng = new Random(42);
            var library = new List<LibrarySpectrum>(precursorCount);
            var windowIds = index.GetWindowIds().OrderBy(w => w).ToList();
            if (windowIds.Count == 0) return library;

            for (int p = 0; p < precursorCount; p++)
            {
                int wid = windowIds[rng.Next(windowIds.Count)];
                var bounds = index.GetWindowBounds(wid);
                double precursorMz = bounds.LowerBound + rng.NextDouble() * (bounds.UpperBound - bounds.LowerBound);

                if (!index.TryGetScanRangeForWindow(wid, out int startScan, out int scanCount) || scanCount == 0)
                    continue;

                float rtStart = index.GetScanRt(startScan);
                float rtEnd = index.GetScanRt(startScan + scanCount - 1);
                double rt = rtStart + rng.NextDouble() * (rtEnd - rtStart);

                int sampleScan = startScan + rng.Next(scanCount);
                var scanMz = index.GetScanMzSpan(sampleScan);

                var fragments = new List<MatchedFragmentIon>(fragmentsPerPrecursor);
                var usedIndices = new HashSet<int>();

                for (int f = 0; f < fragmentsPerPrecursor && usedIndices.Count < scanMz.Length; f++)
                {
                    int peakIdx, attempts = 0;
                    do { peakIdx = rng.Next(scanMz.Length); attempts++; }
                    while (usedIndices.Contains(peakIdx) && attempts < 50);
                    if (usedIndices.Contains(peakIdx)) break;
                    usedIndices.Add(peakIdx);

                    double fragMz = scanMz[peakIdx];
                    var product = new Product(ProductType.b, FragmentationTerminus.N,
                        fragMz - 1.007276, f + 1, 0, 0);
                    fragments.Add(new MatchedFragmentIon(product, fragMz,
                        100.0 + rng.NextDouble() * 9900.0, 1));
                }

                if (fragments.Count >= 3)
                {
                    library.Add(new LibrarySpectrum($"SYNPEPTIDE{p}", precursorMz,
                        2 + rng.Next(3), fragments, rt, p % 10 == 0));
                }
            }
            return library;
        }

        private static MsDataScan[] GenerateSyntheticScans(int windowCount, int scansPerWindow, int peaksPerScan)
        {
            var rng = new Random(42);
            int total = windowCount * scansPerWindow;
            var scans = new MsDataScan[total];
            double windowWidth = 25.0, windowStart = 400.0;
            double windowSpacing = (1200.0 - windowStart) / windowCount;
            int residents = Math.Min(80, peaksPerScan / 2);

            double[][] residentMz = new double[windowCount][];
            for (int w = 0; w < windowCount; w++)
            {
                residentMz[w] = new double[residents];
                for (int f = 0; f < residents; f++)
                    residentMz[w][f] = 100.0 + rng.NextDouble() * 1700.0;
                Array.Sort(residentMz[w]);
            }

            int idx = 0;
            for (int cycle = 0; cycle < scansPerWindow; cycle++)
            {
                double cycleRt = cycle * 0.04;
                for (int w = 0; w < windowCount; w++)
                {
                    double isoCenter = windowStart + w * windowSpacing + windowSpacing / 2.0;
                    double[] mz = new double[peaksPerScan], intensities = new double[peaksPerScan];
                    for (int f = 0; f < residents; f++)
                    {
                        mz[f] = residentMz[w][f] + (rng.NextDouble() - 0.5) * 0.002;
                        intensities[f] = 500.0 + rng.NextDouble() * 5000.0;
                    }
                    for (int p = residents; p < peaksPerScan; p++)
                    {
                        mz[p] = 100.0 + rng.NextDouble() * 1900.0;
                        intensities[p] = rng.NextDouble() * 300.0;
                    }
                    Array.Sort(mz, intensities);

                    scans[idx] = new MsDataScan(new MzSpectrum(mz, intensities, false),
                        idx + 1, 2, true, Polarity.Positive, cycleRt + w * 0.001,
                        new MzRange(100, 2000), "FTMS", MZAnalyzerType.Orbitrap,
                        intensities.Sum(), 20.0, null, $"scan={idx + 1}",
                        isolationMZ: isoCenter, isolationWidth: windowWidth,
                        dissociationType: DissociationType.HCD);
                    idx++;
                }
            }
            return scans;
        }
    }
}