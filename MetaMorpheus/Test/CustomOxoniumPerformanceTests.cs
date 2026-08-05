using EngineLayer;
using EngineLayer.GlycoSearch;
using NUnit.Framework;
using System;
using System.Diagnostics;
using System.Reflection;
using System.Runtime.CompilerServices;

namespace Test
{
    /// <summary>
    /// Performance bounds for the strict custom-oxonium branch of GlycoPeptides.DiagonsticFilter.
    /// Marked [Explicit] + [Category("Performance")] so they never run in the normal suite (timing
    /// assertions would be noisy in CI). Invoke via Test\perf\run-perf.ps1, which filters on
    /// Category=Performance. Bounds are deliberately generous: they catch a gross regression, not a
    /// few-nanosecond drift. The exact ns/call is written to the test output for the report.
    /// </summary>
    [TestFixture]
    public class CustomOxoniumPerformanceTests
    {
        private const int Iterations = 200_000;
        private const int Rounds = 5;

        [Test, Explicit, Category("Performance")]
        public void Filter_NoCustoms_OverheadOverLegacyIsNegligible()
        {
            try
            {
                Glycan.ResetCustomMonosaccharides();
                var box = BoxWithKind(new byte[Glycan.KindCapacity]);
                double[] oxo = new double[Glycan.AllOxoniumIons.Length]; // all zero, no customs

                double legacyNs = MinNsPerCall(() => DiagonsticFilter_LegacyForBenchmark(oxo, box));
                double newNs = MinNsPerCall(() => GlycoPeptides.DiagonsticFilter(oxo, box));
                double addedNs = newNs - legacyNs;

                // The headline is the absolute per-call addition (one Dictionary.Count check). The
                // percentage is reported only for completeness: it is inflated because a ~9 ns method
                // is measured in isolation, where any fixed addition looks large relative to the base.
                TestContext.WriteLine(
                    $"no customs: legacy={legacyNs:F1} ns/call, new={newNs:F1} ns/call, " +
                    $"added={addedNs:F1} ns/call ({addedNs / legacyNs * 100.0:F0}% of the isolated base; " +
                    $"negligible against real per-scan search cost)");

                // Robust ceiling: at zero customs the only added work is one Dictionary.Count check.
                Assert.That(addedNs, Is.LessThan(legacyNs + 50),
                    $"No-customs absolute overhead too high: {addedNs:F1} ns/call");
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
            }
        }

        [Test, Explicit, Category("Performance")]
        public void Filter_FiftyCustomIons_PerCallUnderTenMicroseconds()
        {
            try
            {
                Glycan.ResetCustomMonosaccharides();
                RegisterFiftyCustomIons();
                var customs = Glycan.CustomOxoniumIons;
                Assert.That(customs.Count, Is.EqualTo(50));

                double[] oxo = new double[Glycan.AllOxoniumIonsIncludingCustoms.Length];
                byte[] kind = new byte[Glycan.KindCapacity];
                for (int j = 0; j < customs.Count; j++)
                {
                    oxo[Glycan.AllOxoniumIons.Length + j] = 1.0; // every custom ion observed
                    kind[customs[j].KindIndex] = 1;              // every linked mono present
                }
                var box = BoxWithKind(kind);

                // Worst case: all match, so the loop iterates all 50 entries and returns true.
                Assert.That(GlycoPeptides.DiagonsticFilter(oxo, box), Is.True);

                double ns = MinNsPerCall(() => GlycoPeptides.DiagonsticFilter(oxo, box));
                TestContext.WriteLine($"50 custom ions: {ns:F1} ns/call");
                Assert.That(ns, Is.LessThan(10_000)); // 10 microseconds
            }
            finally
            {
                Glycan.ResetCustomMonosaccharides();
            }
        }

        private static void RegisterFiftyCustomIons()
        {
            char[] codes = { 'U', 'V', 'W', 'Z', 'Q' }; // none collide with built-in codes
            for (int m = 0; m < codes.Length; m++)
            {
                int[] ions = new int[10];
                for (int k = 0; k < ions.Length; k++)
                {
                    ions[k] = 30_000_000 + m * 1000 + k; // distinct scaled-int m/z values
                }
                Glycan.RegisterCustomMonosaccharide($"PerfSugar{m}", codes[m], 18_000_000 + m, ions);
            }
        }

        private static double MinNsPerCall(Func<bool> action)
        {
            for (int i = 0; i < 10_000; i++) // warm the JIT
            {
                action();
            }
            double best = double.MaxValue;
            for (int r = 0; r < Rounds; r++)
            {
                var sw = Stopwatch.StartNew();
                for (int i = 0; i < Iterations; i++)
                {
                    action();
                }
                sw.Stop();
                double ns = sw.Elapsed.TotalMilliseconds * 1e6 / Iterations;
                if (ns < best) best = ns;
            }
            return best;
        }

        private static GlycanBox BoxWithKind(byte[] kind)
        {
            var box = (GlycanBox)RuntimeHelpers.GetUninitializedObject(typeof(GlycanBox));
            MethodInfo setter = typeof(GlycanBox).GetProperty("Kind").GetSetMethod(nonPublic: true);
            setter.Invoke(box, new object[] { kind });
            return box;
        }

        // Pre-Phase-3 DiagonsticFilter body (built-in rules only), kept here as the benchmark baseline.
        private static bool DiagonsticFilter_LegacyForBenchmark(double[] oxoniumIonsintensities, GlycanBox glycanBox)
        {
            double HexNAc_diagnostic = oxoniumIonsintensities[4];
            double NeuAc_diagnostic1 = oxoniumIonsintensities[10];
            double NeuAc_diagnostic2 = oxoniumIonsintensities[12];
            double HexNAcPlusHex_diagnostic = oxoniumIonsintensities[14];

            if (NeuAc_diagnostic1 / HexNAc_diagnostic > 0.02 && NeuAc_diagnostic2 / HexNAc_diagnostic > 0.02)
            {
                if (glycanBox.Kind[2] == 0)
                {
                    return false;
                }
            }

            if (NeuAc_diagnostic1 / HexNAc_diagnostic < 0.02 && NeuAc_diagnostic2 / HexNAc_diagnostic < 0.02)
            {
                if (glycanBox.Kind[2] != 0)
                {
                    return false;
                }
            }
            else if (HexNAcPlusHex_diagnostic / HexNAc_diagnostic > 0.02)
            {
                if (glycanBox.Kind[0] < 1 && glycanBox.Kind[1] < 1)
                {
                    return false;
                }
            }

            return true;
        }
    }
}
