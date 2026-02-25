// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;

namespace Development.Dia
{
    /// <summary>
    /// Benchmarks for Phase 7 library bridge components.
    /// </summary>
    public static class LibraryBridgeBenchmark
    {
        /// <summary>
        /// Builds a realistic DiaScanIndex via DiaScanIndexBuilder for benchmarking.
        /// </summary>
        private static DiaScanIndex BuildBenchmarkIndex(int windowCount, int cyclesPerWindow, int peaksPerScan)
        {
            var scans = new List<MsDataScan>();
            int scanNumber = 1;
            float startMz = 400f;
            float windowWidth = (1200f - 400f) / windowCount;

            for (int cycle = 0; cycle < cyclesPerWindow; cycle++)
            {
                for (int w = 0; w < windowCount; w++)
                {
                    double rt = cycle * 0.6 + w * 0.01;
                    double center = startMz + w * windowWidth + windowWidth / 2;

                    double[] mzs = new double[peaksPerScan];
                    double[] intensities = new double[peaksPerScan];
                    for (int p = 0; p < peaksPerScan; p++)
                    {
                        mzs[p] = 100.0 + p * 0.5;
                        intensities[p] = 1000.0;
                    }

                    var spectrum = new MzSpectrum(mzs, intensities, false);
                    scans.Add(new MsDataScan(
                        massSpectrum: spectrum,
                        oneBasedScanNumber: scanNumber++,
                        msnOrder: 2,
                        isCentroid: true,
                        polarity: Polarity.Positive,
                        retentionTime: rt,
                        scanWindowRange: new MzRange(100, 1500),
                        scanFilter: "",
                        mzAnalyzer: MZAnalyzerType.Orbitrap,
                        totalIonCurrent: intensities.Sum(),
                        injectionTime: 50,
                        noiseData: null,
                        nativeId: $"scan={scanNumber - 1}",
                        isolationMZ: center,
                        isolationWidth: windowWidth,
                        dissociationType: DissociationType.HCD
                    ));
                }
            }

            return DiaScanIndexBuilder.Build(scans.ToArray());
        }

        public static void BenchmarkFindWindow(int windowCount = 50, int lookupCount = 1_000_000)
        {
            Console.WriteLine($"=== FindWindowForPrecursorMz Benchmark ===");
            Console.WriteLine($"Windows: {windowCount}, Lookups: {lookupCount:N0}");

            var index = BuildBenchmarkIndex(windowCount, cyclesPerWindow: 100, peaksPerScan: 200);

            var rng = new Random(42);
            float[] precursors = new float[lookupCount];
            for (int i = 0; i < lookupCount; i++)
                precursors[i] = 380f + (float)(rng.NextDouble() * 840);

            // Warmup
            for (int i = 0; i < 1000; i++)
                index.FindWindowForPrecursorMz(precursors[i % lookupCount]);

            var sw = Stopwatch.StartNew();
            int found = 0;
            for (int i = 0; i < lookupCount; i++)
            {
                if (index.FindWindowForPrecursorMz(precursors[i]) >= 0)
                    found++;
            }
            sw.Stop();

            Console.WriteLine($"Time: {sw.ElapsedMilliseconds} ms");
            Console.WriteLine($"Throughput: {lookupCount / sw.Elapsed.TotalSeconds:N0} lookups/sec");
            Console.WriteLine($"Found: {found:N0} / {lookupCount:N0} ({100.0 * found / lookupCount:F1}%)");
            Console.WriteLine();

            index.Dispose();
        }

        public static void BenchmarkQueryGeneration(int precursorCount = 100_000, int fragmentsPerPrecursor = 12)
        {
            Console.WriteLine($"=== DiaLibraryQueryGenerator.Generate() Benchmark ===");
            Console.WriteLine($"Precursors: {precursorCount:N0}, Fragments/precursor: {fragmentsPerPrecursor}");

            var index = BuildBenchmarkIndex(windowCount: 50, cyclesPerWindow: 100, peaksPerScan: 200);

            var rng = new Random(42);
            var precursors = new List<LibraryPrecursorInput>(precursorCount);
            for (int p = 0; p < precursorCount; p++)
            {
                double precursorMz = 400 + rng.NextDouble() * 800;
                float[] fragMzs = new float[fragmentsPerPrecursor];
                float[] fragIntensities = new float[fragmentsPerPrecursor];
                for (int f = 0; f < fragmentsPerPrecursor; f++)
                {
                    fragMzs[f] = 100f + (float)(rng.NextDouble() * 1500);
                    fragIntensities[f] = (float)(rng.NextDouble() * 10000);
                }
                precursors.Add(new LibraryPrecursorInput(
                    $"PEPTIDE{p}", precursorMz, rng.Next(2, 5),
                    rng.NextDouble() * 60, false, fragMzs, fragIntensities));
            }

            var parameters = new DiaSearchParameters { PpmTolerance = 20f, RtToleranceMinutes = 5f, MinFragmentsRequired = 3 };

            // Warmup
            DiaLibraryQueryGenerator.Generate(precursors.GetRange(0, Math.Min(100, precursors.Count)), index, parameters);

            var sw = Stopwatch.StartNew();
            var result = DiaLibraryQueryGenerator.Generate(precursors, index, parameters);
            sw.Stop();

            Console.WriteLine($"Time: {sw.ElapsedMilliseconds} ms");
            Console.WriteLine($"Precursor throughput: {precursorCount / sw.Elapsed.TotalSeconds:N0}/sec");
            Console.WriteLine($"Query throughput: {result.Queries.Length / sw.Elapsed.TotalSeconds:N0}/sec");
            Console.WriteLine($"Queries: {result.Queries.Length:N0}, Groups: {result.PrecursorGroups.Length:N0}");
            Console.WriteLine($"Skipped (no window): {result.SkippedNoWindow:N0}");
            Console.WriteLine();

            index.Dispose();
        }

        public static void RunAll()
        {
            Console.WriteLine("========================================");
            Console.WriteLine("  Library Bridge Benchmarks (Phase 7)");
            Console.WriteLine("========================================\n");

            BenchmarkFindWindow();
            BenchmarkQueryGeneration(10_000, 12);
            BenchmarkQueryGeneration(100_000, 12);

            Console.WriteLine("Done.");
        }
    }
}
