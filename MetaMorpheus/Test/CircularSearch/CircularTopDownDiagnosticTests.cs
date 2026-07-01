using Chemistry;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test.CircularSearch
{
    /// <summary>
    /// Diagnostic test for top-down circular search of AFYRHTQESG against my.mgf.
    ///
    /// The ring input is "AFYRHTQESG" (N=10). CircularProtein will canonicalize it.
    /// Protease: top-down (no cleavage, no linear products).
    ///
    /// This test does NOT assert a specific result — it reports everything the
    /// engine sees so we can diagnose why 0 PSMs are returned.
    ///
    /// It runs two passes against each candidate scan:
    ///   Pass A — default parameters (peak filtering ON, max 200 peaks/window)
    ///   Pass B — filtering OFF (trimMsMsPeaks=false, no intensity ratio cutoff)
    ///
    /// This directly shows whether peak filtering is discarding the internal
    /// fragment ions needed for scoring.
    /// </summary>
    [TestFixture]
    public static class CircularTopDownDiagnosticTests
    {
        private static string TestDataDir =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");

        private static string MgfPath =>
            Path.Combine(TestDataDir, "my.mgf");

        private static string FastaPath =>
            Path.Combine(TestDataDir, "AFYRHTQESG_cyclic.fasta");

        private static readonly double Proton = 1.007276;

        private static void EnsureFastaExists()
        {
            if (File.Exists(FastaPath)) return;
            File.WriteAllText(FastaPath,
                ">sp|CYC0003R00|CYC0003R00_SYNTH Cyclic rotation 0 of AFYRHTQESG " +
                "OS=Synthetic construct OX=32630 GN=CYC0003 PE=4 SV=1\n" +
                "AFYRHTQESG\n");
        }

        /// <summary>
        /// Builds CommonParameters for one diagnostic pass.
        /// </summary>
        private static CommonParameters MakeParams(bool filterPeaks) =>
            new CommonParameters(
                dissociationType: DissociationType.HCD,
                precursorMassTolerance: new PpmTolerance(20),
                productMassTolerance: new PpmTolerance(20),
                scoreCutoff: 1,
                numberOfPeaksToKeepPerWindow: filterPeaks ? 200 : int.MaxValue,
                minimumAllowedIntensityRatioToBasePeak: filterPeaks ? 0.01 : 0.0,
                trimMs1Peaks: false,
                trimMsMsPeaks: !filterPeaks ? false : true,
                digestionParams: new DigestionParams(
                    protease: "top-down",
                    maxMissedCleavages: 0,
                    minPeptideLength: 1));

        [Test]
        public static void Diagnostic_TopDown_AFYRHTQESG_InspectInternalIonMatching()
        {
            Assert.That(File.Exists(MgfPath),
                $"MGF not found. Place 'my.mgf' in {TestDataDir}.");

            EnsureFastaExists();

            // ── Step 1: Canonicalize and digest ──────────────────────────────
            var protein = new CircularProtein("AFYRHTQESG", "CYC0003R00");
            Console.WriteLine("=== PROTEIN ===");
            Console.WriteLine($"Input sequence:     AFYRHTQESG");
            Console.WriteLine($"Canonical sequence: {protein.BaseSequence}");
            Console.WriteLine($"Ring length N:      {protein.BaseSequence.Length}");

            var commonParams = MakeParams(filterPeaks: false);

            var digestProducts = protein
                .Digest(commonParams.DigestionParams,
                    new List<Modification>(),
                    new List<Modification>())
                .ToList();

            Console.WriteLine($"\n=== DIGESTION (top-down) ===");
            Console.WriteLine($"Total products: {digestProducts.Count}");

            var circularPeptide = digestProducts
                .OfType<CircularPeptideWithSetModifications>()
                .SingleOrDefault();

            if (circularPeptide == null)
            {
                Console.WriteLine("ERROR: No CircularPeptideWithSetModifications produced.");
                Assert.Fail("Top-down digestion produced no circular product.");
                return;
            }

            Console.WriteLine($"Circular peptide:   {circularPeptide.BaseSequence}");
            Console.WriteLine($"MonoisotopicMass:   {circularPeptide.MonoisotopicMass:F6} Da");
            Console.WriteLine($"Precursor m/z (+1): {circularPeptide.MonoisotopicMass + Proton:F6}");
            Console.WriteLine($"Precursor m/z (+2): {(circularPeptide.MonoisotopicMass + 2 * Proton) / 2:F6}");
            Console.WriteLine($"Precursor m/z (+3): {(circularPeptide.MonoisotopicMass + 3 * Proton) / 3:F6}");

            // ── Step 2: Generate all internal fragments ───────────────────────
            var fragmentProducts = new List<Product>();
            circularPeptide.FragmentInternally(
                DissociationType.HCD,
                minLengthOfFragments: 2,
                fragmentProducts);

            Console.WriteLine($"\n=== THEORETICAL INTERNAL FRAGMENTS (minLength=2, HCD) ===");
            Console.WriteLine($"Total: {fragmentProducts.Count}");
            Console.WriteLine($"{"Fragment",-20} {"NeutralMass":>14} {"m/z (+1)":>12}");
            Console.WriteLine(new string('-', 50));
            foreach (var frag in fragmentProducts.OrderBy(f => f.NeutralMass))
            {
                string label = $"[{frag.FragmentNumber}-{frag.SecondaryFragmentNumber}]";
                Console.WriteLine(
                    $"{label,-20} {frag.NeutralMass,14:F6} {frag.NeutralMass + Proton,12:F6}");
            }

            // ── Step 3: Load scans — two versions, filtered and unfiltered ────
            Console.WriteLine("\n=== MGF SCANS ===");

            var paramsFiltered = MakeParams(filterPeaks: true);
            var paramsUnfiltered = MakeParams(filterPeaks: false);

            var msDataFileFiltered = new Readers.Mgf(MgfPath).LoadAllStaticData();
            var msDataFileUnfiltered = new Readers.Mgf(MgfPath).LoadAllStaticData();

            var scansFiltered = MetaMorpheusTask
                .GetMs2Scans(msDataFileFiltered, MgfPath, paramsFiltered)
                .OrderBy(s => s.PrecursorMass).ToArray();
            var scansUnfiltered = MetaMorpheusTask
                .GetMs2Scans(msDataFileUnfiltered, MgfPath, paramsUnfiltered)
                .OrderBy(s => s.PrecursorMass).ToArray();

            Console.WriteLine($"Total MS2 scans loaded: {scansFiltered.Length}");

            if (scansFiltered.Length == 0)
            {
                Console.WriteLine("ERROR: No MS2 scans found in MGF.");
                Assert.Fail("No MS2 scans in my.mgf.");
                return;
            }

            // ── Step 4: Precursor mass check ──────────────────────────────────
            Console.WriteLine($"\n=== PRECURSOR MASS CHECK (±20 ppm) ===");
            double theorMass = circularPeptide.MonoisotopicMass;
            double tolDa = theorMass * 20.0 / 1e6;
            Console.WriteLine($"Theoretical mass:  {theorMass:F6} Da");
            Console.WriteLine($"Tolerance window:  [{theorMass - tolDa:F6}, {theorMass + tolDa:F6}] Da");
            Console.WriteLine();
            Console.WriteLine(
                $"{"Scan",-8} {"PrecursorMass":>14} {"Error(ppm)":>12} " +
                $"{"InWindow?":>10} {"ExpFrags(filt)":>14} {"ExpFrags(unfilt)":>16}");
            Console.WriteLine(new string('-', 80));

            var matchingScansFiltered = new List<Ms2ScanWithSpecificMass>();
            var matchingScansUnfiltered = new List<Ms2ScanWithSpecificMass>();

            for (int i = 0; i < scansFiltered.Length; i++)
            {
                var sf = scansFiltered[i];
                var su = scansUnfiltered[i];
                double errorPpm = (sf.PrecursorMass - theorMass) / theorMass * 1e6;
                bool inWindow = Math.Abs(errorPpm) <= 20.0;
                int expFiltCount = sf.ExperimentalFragments?.Length ?? 0;
                int expUnfiltCount = su.ExperimentalFragments?.Length ?? 0;

                Console.WriteLine(
                    $"{sf.OneBasedScanNumber,-8} {sf.PrecursorMass,14:F6} " +
                    $"{errorPpm,12:F2} {(inWindow ? "YES" : "no"),10} " +
                    $"{expFiltCount,14} {expUnfiltCount,16}");

                if (inWindow)
                {
                    matchingScansFiltered.Add(sf);
                    matchingScansUnfiltered.Add(su);
                }
            }

            Console.WriteLine($"\nScans within precursor tolerance: {matchingScansFiltered.Count}");

            if (matchingScansFiltered.Count == 0)
            {
                Console.WriteLine(
                    "\nDIAGNOSIS: No scan within 20 ppm of theoretical mass." +
                    "\nPossible causes:" +
                    "\n  1. Different charge state — check m/z vs neutral mass encoding." +
                    "\n  2. Post-translational or cyclization modification on the ring." +
                    "\n  3. MGF PEPMASS is m/z and charge is unknown/wrong." +
                    "\n  4. Precursor tolerance needs widening.");

                var closest = scansFiltered
                    .OrderBy(s => Math.Abs(s.PrecursorMass - theorMass))
                    .First();
                double closestPpm =
                    (closest.PrecursorMass - theorMass) / theorMass * 1e6;
                Console.WriteLine(
                    $"\nClosest scan: #{closest.OneBasedScanNumber}  " +
                    $"precursorMass={closest.PrecursorMass:F6}  " +
                    $"error={closestPpm:F2} ppm");

                Assert.Pass("Diagnostic complete — precursor mass mismatch. See console.");
                return;
            }

            // ── Step 5: Fragment matching — filtered vs unfiltered ────────────
            Console.WriteLine(
                $"\n=== FRAGMENT MATCHING: FILTERED vs UNFILTERED ===");
            Console.WriteLine(
                $"{"Scan",-8} {"Matched(filt)":>14} {"Score(filt)":>12} " +
                $"{"Matched(unfilt)":>16} {"Score(unfilt)":>14}");
            Console.WriteLine(new string('-', 70));

            int bestCountUnfilt = 0;
            Ms2ScanWithSpecificMass bestScanUnfilt = null;
            List<MatchedFragmentIon> bestIonsUnfilt = null;

            for (int i = 0; i < matchingScansFiltered.Count; i++)
            {
                var sf = matchingScansFiltered[i];
                var su = matchingScansUnfiltered[i];

                var matchedFilt = MetaMorpheusEngine.MatchFragmentIons(
                    sf, fragmentProducts, paramsFiltered);
                var matchedUnfilt = MetaMorpheusEngine.MatchFragmentIons(
                    su, fragmentProducts, paramsUnfiltered);

                double scoreFilt = MetaMorpheusEngine.CalculatePeptideScore(
                    sf.TheScan, matchedFilt);
                double scoreUnfilt = MetaMorpheusEngine.CalculatePeptideScore(
                    su.TheScan, matchedUnfilt);

                Console.WriteLine(
                    $"{sf.OneBasedScanNumber,-8} " +
                    $"{matchedFilt.Count,14} {scoreFilt,12:F3} " +
                    $"{matchedUnfilt.Count,16} {scoreUnfilt,14:F3}");

                if (matchedUnfilt.Count > bestCountUnfilt)
                {
                    bestCountUnfilt = matchedUnfilt.Count;
                    bestScanUnfilt = su;
                    bestIonsUnfilt = matchedUnfilt;
                }
            }

            // ── Step 6: Detailed report for best unfiltered scan ──────────────
            if (bestScanUnfilt != null)
            {
                Console.WriteLine(
                    $"\n=== BEST SCAN (unfiltered): #{bestScanUnfilt.OneBasedScanNumber} ===");
                Console.WriteLine($"Matched ions: {bestIonsUnfilt.Count}");

                if (bestIonsUnfilt.Count > 0)
                {
                    Console.WriteLine(
                        $"\n{"Fragment",-20} {"TheoMass":>12} {"ObsMass":>12} " +
                        $"{"Error(ppm)":>12} {"Intensity":>12}");
                    Console.WriteLine(new string('-', 72));

                    foreach (var ion in bestIonsUnfilt.OrderBy(
                        i => i.NeutralTheoreticalProduct.FragmentNumber))
                    {
                        string label =
                            $"[{ion.NeutralTheoreticalProduct.FragmentNumber}-" +
                            $"{ion.NeutralTheoreticalProduct.SecondaryFragmentNumber}]";
                        double theoMz = ion.NeutralTheoreticalProduct.NeutralMass;
                        double obsMass = ion.Mz.ToMass(ion.Charge);
                        double errPpm = (obsMass - theoMz) / theoMz * 1e6;
                        Console.WriteLine(
                            $"{label,-20} {theoMz,12:F6} {obsMass,12:F6} " +
                            $"{errPpm,12:F2} {ion.Intensity,12:F0}");
                    }
                }
                else
                {
                    // Report what IS in the spectrum vs what we expect
                    Console.WriteLine(
                        "\nDIAGNOSIS: Precursor matches but no fragment ions match." +
                        "\nComparing theoretical vs experimental masses:");

                    Console.WriteLine(
                        $"\nTheoretical internal fragment m/z values (charge +1):");
                    foreach (var f in fragmentProducts.OrderBy(f => f.NeutralMass))
                        Console.WriteLine(
                            $"  [{f.FragmentNumber}-{f.SecondaryFragmentNumber}]  " +
                            $"{f.NeutralMass + Proton:F4}");

                    Console.WriteLine(
                        $"\nExperimental fragment neutral masses (scan " +
                        $"#{bestScanUnfilt.OneBasedScanNumber}, unfiltered):");
                    if (bestScanUnfilt.ExperimentalFragments?.Length > 0)
                    {
                        foreach (var ef in bestScanUnfilt.ExperimentalFragments
                            .OrderBy(e => e.MonoisotopicMass))
                            Console.WriteLine(
                                $"  {ef.MonoisotopicMass:F4} Da  " +
                                $"charge={ef.Charge}  " +
                                $"intensity={ef.TotalIntensity:F0}");
                    }
                    else
                    {
                        Console.WriteLine(
                            "  (none — deconvolution produced no fragments even " +
                            "without filtering)");
                    }

                    Console.WriteLine(
                        "\nPossible causes:" +
                        "\n  1. Widen product mass tolerance (try 50 or 100 ppm)." +
                        "\n  2. Internal fragment masses require a different dissociation " +
                            "type (try CID)." +
                        "\n  3. The experimental peaks are present but at charge > 1 — " +
                            "deconvolution may be needed.");
                }
            }

            Assert.Pass("Diagnostic complete — see console output for full analysis.");
        }
    }
}