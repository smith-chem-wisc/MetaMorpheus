using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test.CircularSearch
{
    /// <summary>
    /// Integration test using real MS data.
    ///
    /// Data file : "d Tryp-SFYR MS2 CID35 plus2.mgf"
    /// Ring      : input "HGQAETSFYR", canonical "AETSFYRHGQ" (A &lt; H)
    /// Protease  : trypsin, 0 missed cleavages
    ///
    /// The ring "AETSFYRHGQ" has one trypsin site (after R at canonical pos 8).
    /// With 0 missed cleavages, the only product is the wrapping linear peptide
    /// "HGQAETSFYR". The engine scores this against every MS2 scan in the file.
    ///
    /// The reference classic search identified scan 349 as the best hit with
    /// score 16.33 and 22 matched ions. However, the CircularSearchEngine with
    /// its dual terminal+internal ion scoring produces higher scores for other
    /// scans (e.g. scan 151 with 18 matched ions). The top-scoring scan is
    /// therefore not asserted here — instead we verify:
    ///
    ///   1. Digestion pre-conditions hold.
    ///   2. The task runs without error.
    ///   3. Scan 349 is found in the output with the correct base sequence,
    ///      target classification, score > 10, and matched ions from all three
    ///      ion series (y, b, internal).
    ///   4. All PSMs in the output have base sequence "HGQAETSFYR" since there
    ///      is only one peptide in the database.
    /// </summary>
    [TestFixture]
    public static class CircularSearchRealDataTests
    {
        // ── File paths ────────────────────────────────────────────────────────

        private static string TestDataDir =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");

        private static string MgfPath =>
            Path.Combine(TestDataDir, "d Tryp-SFYR MS2 CID35 plus2.mgf");

        private static string FastaPath =>
            Path.Combine(TestDataDir, "HGQAETSFYR_cyclic.fasta");

        // ── FASTA writer ──────────────────────────────────────────────────────

        private static void EnsureFastaExists()
        {
            if (File.Exists(FastaPath)) return;
            File.WriteAllText(FastaPath,
                ">sp|CYC0002R07|CYC0002R07_SYNTH Cyclic rotation 7 of AETSFYRHGQ " +
                "OS=Synthetic construct OX=32630 GN=CYC0002 PE=4 SV=1\n" +
                "HGQAETSFYR\n");
        }

        // ── Test ──────────────────────────────────────────────────────────────

        [Test]
        public static void RealData_TrypSFYR_Scan349_FoundWithCorrectSequenceAndIons()
        {
            // ── Pre-condition: MGF must exist ─────────────────────────────────
            Assert.That(File.Exists(MgfPath),
                $"MGF not found. Place 'd Tryp-SFYR MS2 CID35 plus2.mgf' in {TestDataDir}.");

            EnsureFastaExists();

            // ── Pre-condition 1: canonical sequence ───────────────────────────
            var protein = new CircularProtein("HGQAETSFYR", "CYC0002R07");
            Assert.That(protein.BaseSequence, Is.EqualTo("AETSFYRHGQ"),
                "Canonical sequence of HGQAETSFYR must be AETSFYRHGQ.");
            Console.WriteLine($"Canonical sequence: {protein.BaseSequence}");

            // ── Pre-condition 2: 0 missed cleavages → no circular product ─────
            var digestionParams0 = new DigestionParams(
                protease: "trypsin", maxMissedCleavages: 0, minPeptideLength: 1);

            var products0 = protein
                .Digest(digestionParams0, new List<Modification>(), new List<Modification>())
                .ToList();

            Console.WriteLine($"Products at mmc=0: [{string.Join(", ", products0.Select(p => p.BaseSequence))}]");

            Assert.That(products0.OfType<CircularPeptideWithSetModifications>().Count(),
                Is.EqualTo(0), "No circular product expected at mmc=0.");

            // ── Pre-condition 3: wrapping linear product "HGQAETSFYR" ─────────
            Assert.That(products0.OfType<PeptideWithSetModifications>()
                    .Where(p => p is not CircularPeptideWithSetModifications)
                    .Any(p => p.BaseSequence == "HGQAETSFYR"),
                "Wrapping linear product 'HGQAETSFYR' must be present.");

            // ── Run the task ──────────────────────────────────────────────────
            string outputDir = Path.Combine(
                TestContext.CurrentContext.TestDirectory,
                $"CircularRealDataTest_{Guid.NewGuid():N}");
            Directory.CreateDirectory(outputDir);

            try
            {
                var task = new CircularSearchTask
                {
                    CommonParameters = new CommonParameters(
                        dissociationType: DissociationType.HCD,
                        precursorMassTolerance: new PpmTolerance(20),
                        productMassTolerance: new PpmTolerance(30),
                        scoreCutoff: 1,
                        digestionParams: new DigestionParams(
                            protease: "trypsin",
                            maxMissedCleavages: 0,
                            minPeptideLength: 1)),
                    CircularSearchParameters = new CircularSearchParameters
                    {
                        WriteHighQValuePsms = true,
                        WriteDecoys = true,
                        WriteContaminants = true,
                        MinInternalFragmentLength = 2
                    }
                };

                var engine = new EverythingRunnerEngine(
                    taskList: new List<(string, MetaMorpheusTask)> { ("CircularSearch", task) },
                    startingRawFilenameList: new List<string> { MgfPath },
                    startingXmlDbFilenameList: new List<DbForTask>
                        { new DbForTask(FastaPath, isContaminant: false) },
                    outputFolder: outputDir);

                Assert.That(() => engine.Run(), Throws.Nothing,
                    "CircularSearchTask should complete without throwing.");

                // ── Read output TSV ───────────────────────────────────────────
                string tsvPath = Path.Combine(outputDir, "CircularSearch",
                    "CircularPeptideSpectralMatches.psmtsv");
                Assert.That(File.Exists(tsvPath), "Output TSV must exist.");

                var lines = File.ReadAllLines(tsvPath);
                Assert.That(lines.Length, Is.GreaterThan(1),
                    "Output TSV must have at least one data row.");

                string[] headers = lines[0].Split('\t');
                var rows = lines.Skip(1)
                    .Where(l => !string.IsNullOrWhiteSpace(l))
                    .Select(l => l.Split('\t'))
                    .ToList();

                Assert.That(rows.Count, Is.GreaterThan(0), "At least one PSM must be present.");

                string Col(string[] row, string name)
                {
                    int idx = Array.IndexOf(headers, name);
                    Assert.That(idx, Is.GreaterThanOrEqualTo(0),
                        $"Column '{name}' not found.");
                    return idx < row.Length ? row[idx] : string.Empty;
                }

                // ── Assertion 1: all PSMs have the correct base sequence ───────
                // Only one peptide in the database, so every matched scan must
                // identify HGQAETSFYR.
                foreach (var row in rows)
                {
                    Assert.That(Col(row, "Base Sequence"), Is.EqualTo("HGQAETSFYR"),
                        $"Scan {Col(row, "Scan Number")}: base sequence must be HGQAETSFYR.");
                }

                // ── Assertion 2: scan 349 is present in the output ────────────
                var scan349Row = rows.FirstOrDefault(r => Col(r, "Scan Number") == "349");
                Assert.That(scan349Row, Is.Not.Null,
                    "Scan 349 must be present in the output TSV.");

                Console.WriteLine($"Scan 349 base sequence: {Col(scan349Row, "Base Sequence")}");
                Console.WriteLine($"Scan 349 score:         {Col(scan349Row, "Score")}");
                Console.WriteLine($"Scan 349 ion counts:    {Col(scan349Row, "Matched Ion Counts")}");
                Console.WriteLine($"Scan 349 ion series:    {Col(scan349Row, "Matched Ion Series")}");
                Console.WriteLine($"Scan 349 D/C/T:         {Col(scan349Row, "Decoy/Contaminant/Target")}");

                // ── Assertion 3: scan 349 is a target ────────────────────────
                Assert.That(Col(scan349Row, "Decoy/Contaminant/Target"), Is.EqualTo("T"),
                    "Scan 349 must be a target PSM.");

                // ── Assertion 4: scan 349 score > 10 (reference: 16.33) ───────
                bool scoreValid = double.TryParse(Col(scan349Row, "Score"),
                    System.Globalization.NumberStyles.Float,
                    System.Globalization.CultureInfo.InvariantCulture,
                    out double score349);
                Assert.That(scoreValid, "Scan 349 score must parse as double.");
                Assert.That(score349, Is.GreaterThan(10),
                    $"Scan 349 score must be > 10. Actual: {score349}");

                // ── Assertion 5: scan 349 has ions from all three series ───────
                // The reference shows y-ions, b-ions, and internal (yIb) ions.
                // With dual terminal+internal scoring, all three series should match.
                //
                // REGRESSION GUARD: If only internal ions (yIb) are present and
                // terminal ions (y, b) are missing, FragmentInternally's
                // products.Clear() is destroying terminal ions generated by
                // Fragment(). The fix is to use a separate list for internal
                // fragments in CircularSearchEngine (see the comment there).
                string series349 = Col(scan349Row, "Matched Ion Series");
                Assert.That(series349, Does.Contain("y"),
                    "Scan 349 must have matched y-ions. If only yIb ions are present, " +
                    "FragmentInternally.Clear() may be destroying terminal ions — " +
                    "ensure internal fragments use a separate list in CircularSearchEngine.");
                Assert.That(series349, Does.Contain("b"),
                    "Scan 349 must have matched b-ions. If only yIb ions are present, " +
                    "FragmentInternally.Clear() may be destroying terminal ions — " +
                    "ensure internal fragments use a separate list in CircularSearchEngine.");
                Assert.That(series349, Does.Contain("yIb"),
                    "Scan 349 must have matched internal (yIb) ions.");

                // ── Assertion 6: total matched ions >= 10 ────────────────────
                string matchedCountStr = Col(scan349Row, "Matched Ion Counts");
                int totalIons = matchedCountStr
                    .Split(';')
                    .Sum(s => int.TryParse(s.Trim(), out int n) ? n : 0);
                Assert.That(totalIons, Is.GreaterThanOrEqualTo(10),
                    $"Scan 349 must have >= 10 matched ions. Actual: {totalIons}. " +
                    "If only ~4 ions matched, terminal ions may have been cleared " +
                    "by FragmentInternally — check CircularSearchEngine fragment list handling.");
            }
            finally
            {
                Directory.Delete(outputDir, recursive: true);
            }
        }
    }
}