// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace Development.Dia
{
    public static class Program
    {
        public static void Main(string[] args)
        {
            Console.WriteLine("================================================================");
            Console.WriteLine("  MetaMorpheus DIA Engine — Development Benchmarks");
            Console.WriteLine("================================================================");
            Console.WriteLine($"  Time: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
            Console.WriteLine($"  CLR:  {System.Runtime.InteropServices.RuntimeInformation.FrameworkDescription}");
            Console.WriteLine();

            bool runAll = args.Length == 0;
            bool runEngine = runAll || Array.Exists(args, a => a.Equals("--engine", StringComparison.OrdinalIgnoreCase));
            bool runBridge = runAll || Array.Exists(args, a => a.Equals("--bridge", StringComparison.OrdinalIgnoreCase));

            if (args.Length > 0 && Array.Exists(args, a => a.Equals("--help", StringComparison.OrdinalIgnoreCase)))
            {
                Console.WriteLine("Usage: DiaDevBenchmarks [--engine] [--bridge] [--help]");
                Console.WriteLine("  No options = run all benchmarks.");
                return;
            }

            if (runEngine) DiaEngineBenchmark.RunAll();
            if (runBridge) LibraryBridgeBenchmark.RunAll();
            DiaEndToEndBenchmark.RunAll();
            // Real-data end-to-end benchmark with .msp libraries + mzML
            DiaRealDataBenchmark.RunAll();
            Console.WriteLine("All benchmarks complete.");
        }
    }
}
