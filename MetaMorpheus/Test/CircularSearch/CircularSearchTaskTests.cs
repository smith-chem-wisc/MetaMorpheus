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
    /// End-to-end tests for <see cref="CircularSearchTask"/> using synthetic .mgf
    /// and .fasta files written to a temp directory.
    ///
    /// Each test writes its own files, runs the task, and inspects the output TSV.
    /// All temp files are cleaned up in a finally block.
    /// </summary>
    [TestFixture]
    public static class CircularSearchTaskTests
    {
        // ── Constants ─────────────────────────────────────────────────────────

        private const string OutputFileName = "CircularPeptideSpectralMatches.psmtsv";
        private static readonly double Proton = 1.007276;

        // ── Infrastructure helpers ────────────────────────────────────────────

        /// <summary>
        /// Creates a unique temp directory for one test run and returns its path.
        /// </summary>
        private static string MakeTempDir(string testName)
        {
            string dir = Path.Combine(
                TestContext.CurrentContext.TestDirectory,
                $"CircularTaskTest_{testName}_{Guid.NewGuid():N}");
            Directory.CreateDirectory(dir);
            return dir;
        }

        /// <summary>
        /// Writes a UniProt-style FASTA file recognised by ProteinDbLoader.
        /// Format: >sp|ACCESSION|NAME Description OS=Synthetic construct OX=32630 GN=ACC PE=4 SV=1
        /// </summary>
        private static string WriteFasta(
            string path,
            IEnumerable<(string accession, string name, string sequence)> entries)
        {
            using var sw = new StreamWriter(path);
            foreach (var (acc, name, seq) in entries)
            {
                sw.WriteLine($">sp|{acc}|{name} Synthetic cyclic peptide OS=Synthetic construct OX=32630 GN={acc} PE=4 SV=1");
                sw.WriteLine(seq);
            }
            return path;
        }

        /// <summary>
        /// Writes a single-scan MGF whose precursor and fragments exactly match
        /// the circular peptide produced by digesting <paramref name="canonicalSequence"/>
        /// with trypsin at the given maxMissedCleavages.
        /// </summary>
        private static string WriteMgfForSequence(
            string mgfPath,
            string canonicalSequence,
            int maxMissedCleavages,
            out double precursorMass)
        {
            var protein = new CircularProtein(canonicalSequence, "tmp");
            var digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: maxMissedCleavages,
                minPeptideLength: 1);

            var circularPeptide = protein
                .Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .OfType<CircularPeptideWithSetModifications>()
                .FirstOrDefault()
                ?? throw new InvalidOperationException(
                    $"No CircularPeptideWithSetModifications for '{canonicalSequence}' " +
                    $"with maxMissedCleavages={maxMissedCleavages}.");

            precursorMass = circularPeptide.MonoisotopicMass;

            var fragmentProducts = new List<Omics.Fragmentation.Product>();
            circularPeptide.FragmentInternally(
                DissociationType.HCD,
                minLengthOfFragments: 3,
                fragmentProducts);

            using var sw = new StreamWriter(mgfPath);
            sw.WriteLine("BEGIN IONS");
            sw.WriteLine(FormattableString.Invariant($"PEPMASS={precursorMass + Proton:F6}"));
            sw.WriteLine("CHARGE=1+");
            sw.WriteLine("SCANS=1");
            sw.WriteLine("RTINSECONDS=60");
            foreach (double nm in fragmentProducts.Select(p => p.NeutralMass).Distinct())
                sw.WriteLine(FormattableString.Invariant($"{nm + Proton:F6} 1000000"));
            sw.WriteLine("END IONS");

            return mgfPath;
        }

        /// <summary>
        /// Constructs and runs an <see cref="EverythingRunnerEngine"/> for a single
        /// <see cref="CircularSearchTask"/> against one MGF and one or more FASTA files.
        /// Output lands in <paramref name="outputDir"/>.
        /// </summary>
        private static void RunTask(
            string outputDir,
            List<DbForTask> databases,
            string mgfPath,
            CircularSearchTask task)
        {
            var engine = new EverythingRunnerEngine(
                taskList: new List<(string, MetaMorpheusTask)> { ("CircularSearch", task) },
                startingRawFilenameList: new List<string> { mgfPath },
                startingXmlDbFilenameList: databases,
                outputFolder: outputDir);

            engine.Run();
        }

        /// <summary>
        /// Builds a default <see cref="CircularSearchTask"/> with scoreCutoff=1
        /// so synthetic scans with few matched ions still produce PSMs.
        /// </summary>
        private static CircularSearchTask MakeTask(
            CircularSearchParameters searchParams = null) =>
            new CircularSearchTask
            {
                CommonParameters = new CommonParameters(
                    dissociationType: DissociationType.HCD,
                    precursorMassTolerance: new PpmTolerance(5),
                    productMassTolerance: new PpmTolerance(20),
                    scoreCutoff: 1,
                    digestionParams: new DigestionParams(
                        protease: "trypsin",
                        maxMissedCleavages: 2,
                        minPeptideLength: 1)),
                CircularSearchParameters = searchParams ?? new CircularSearchParameters()
            };

        /// <summary>
        /// Reads the output TSV and returns header columns and data rows.
        /// Asserts the file exists before reading.
        /// </summary>
        private static (string[] headers, List<string[]> rows) ReadOutputTsv(string outputDir)
        {
            // EverythingRunnerEngine places task output in a subfolder named after the task.
            string taskSubfolder = Path.Combine(outputDir, "CircularSearch");
            string path = Path.Combine(taskSubfolder, OutputFileName);

            Assert.That(File.Exists(path),
                $"Expected output file not found: {path}");

            var lines = File.ReadAllLines(path);
            Assert.That(lines.Length, Is.GreaterThan(0),
                "Output TSV must have at least a header row.");

            string[] headers = lines[0].Split('\t');
            var rows = lines.Skip(1)
                .Where(l => !string.IsNullOrWhiteSpace(l))
                .Select(l => l.Split('\t'))
                .ToList();

            return (headers, rows);
        }

        /// <summary>
        /// Returns the value of a named column for a given data row.
        /// </summary>
        private static string GetColumn(string[] headers, string[] row, string columnName)
        {
            int idx = Array.IndexOf(headers, columnName);
            Assert.That(idx, Is.GreaterThanOrEqualTo(0),
                $"Column '{columnName}' not found in TSV headers. " +
                $"Available: [{string.Join(", ", headers)}]");
            return idx < row.Length ? row[idx] : string.Empty;
        }

        // ── Tests ─────────────────────────────────────────────────────────────

        /// <summary>
        /// Verifies that the task runs end-to-end without error and produces
        /// CircularPeptideSpectralMatches.psmtsv in the output folder.
        ///
        /// Setup: one synthetic FASTA entry ("AFYRHTQESG"), one MGF scan built from
        /// the exact fragments of the circular peptide produced by trypsin digestion.
        /// </summary>
        [Test]
        public static void TaskRuns_ProducesOutputFile()
        {
            string dir = MakeTempDir(nameof(TaskRuns_ProducesOutputFile));
            try
            {
                var protein = new CircularProtein("AFYRHTQESG", "CYC0003R00");
                string canonical = protein.BaseSequence;

                string mgfPath = WriteMgfForSequence(
                    Path.Combine(dir, "synthetic.mgf"), canonical,
                    maxMissedCleavages: 2, out _);

                string fastaPath = WriteFasta(
                    Path.Combine(dir, "synthetic.fasta"),
                    new[] { ("CYC0003R00", "CYC0003R00_SYNTH", canonical) });

                Assert.That(
                    () => RunTask(dir,
                        new List<DbForTask> { new DbForTask(fastaPath, isContaminant: false) },
                        mgfPath,
                        MakeTask()),
                    Throws.Nothing,
                    "CircularSearchTask should complete without throwing.");

                string taskSubfolder = Path.Combine(dir, "CircularSearch");
                string outputPath = Path.Combine(taskSubfolder, OutputFileName);
                Assert.That(File.Exists(outputPath),
                    $"Output file '{OutputFileName}' must exist after task completes.");
            }
            finally { Directory.Delete(dir, recursive: true); }
        }

        /// <summary>
        /// Verifies that WriteDecoys=false produces a TSV with no decoy rows.
        ///
        /// The "Decoy/Contaminant/Target" column contains "D" for decoys.
        /// With WriteDecoys=false, no such rows should appear.
        /// </summary>
        [Test]
        public static void WriteDecoysFalse_NoDecoyRowsInOutput()
        {
            string dir = MakeTempDir(nameof(WriteDecoysFalse_NoDecoyRowsInOutput));
            try
            {
                var protein = new CircularProtein("AFYRHTQESG", "CYC0003R00");
                string canonical = protein.BaseSequence;

                string mgfPath = WriteMgfForSequence(
                    Path.Combine(dir, "synthetic.mgf"), canonical,
                    maxMissedCleavages: 2, out _);

                string fastaPath = WriteFasta(
                    Path.Combine(dir, "synthetic.fasta"),
                    new[] { ("CYC0003R00", "CYC0003R00_SYNTH", canonical) });

                var searchParams = new CircularSearchParameters { WriteDecoys = false };

                RunTask(dir,
                    new List<DbForTask> { new DbForTask(fastaPath, isContaminant: false) },
                    mgfPath,
                    MakeTask(searchParams));

                var (headers, rows) = ReadOutputTsv(dir);

                foreach (var row in rows)
                {
                    string value = GetColumn(headers, row, "Decoy/Contaminant/Target");
                    Assert.That(value, Is.Not.EqualTo("D"),
                        "No decoy rows should appear when WriteDecoys=false.");
                }
            }
            finally { Directory.Delete(dir, recursive: true); }
        }

        /// <summary>
        /// Verifies that WriteContaminants=false produces a TSV with no contaminant rows.
        ///
        /// A contaminant FASTA is loaded via a DbForTask with isContaminant=true.
        /// With WriteContaminants=false, no rows with "C" in the
        /// Decoy/Contaminant/Target column should appear.
        /// </summary>
        [Test]
        public static void WriteContaminantsFalse_NoContaminantRowsInOutput()
        {
            string dir = MakeTempDir(nameof(WriteContaminantsFalse_NoContaminantRowsInOutput));
            try
            {
                var protein = new CircularProtein("AFYRHTQESG", "CYC0003R00");
                string canonical = protein.BaseSequence;

                string mgfPath = WriteMgfForSequence(
                    Path.Combine(dir, "synthetic.mgf"), canonical,
                    maxMissedCleavages: 2, out _);

                // Target FASTA
                string targetFasta = WriteFasta(
                    Path.Combine(dir, "target.fasta"),
                    new[] { ("CYC0003R00", "CYC0003R00_SYNTH", canonical) });

                // Contaminant FASTA — same sequence, different accession
                string contamFasta = WriteFasta(
                    Path.Combine(dir, "contaminant.fasta"),
                    new[] { ("CONTAM001", "CONTAM001_SYNTH", canonical) });

                var searchParams = new CircularSearchParameters { WriteContaminants = false };

                RunTask(dir,
                    new List<DbForTask>
                    {
                        new DbForTask(targetFasta, isContaminant: false),
                        new DbForTask(contamFasta,  isContaminant: true)
                    },
                    mgfPath,
                    MakeTask(searchParams));

                var (headers, rows) = ReadOutputTsv(dir);

                foreach (var row in rows)
                {
                    string value = GetColumn(headers, row, "Decoy/Contaminant/Target");
                    Assert.That(value, Is.Not.EqualTo("C"),
                        "No contaminant rows should appear when WriteContaminants=false.");
                }
            }
            finally { Directory.Delete(dir, recursive: true); }
        }

        /// <summary>
        /// Verifies that WriteHighQValuePsms=false produces a TSV where every row
        /// has QValue ≤ 0.01.
        ///
        /// With a single-scan synthetic file there is no FDR pressure, so the one
        /// matching PSM should have QValue=0. The test verifies that no row has
        /// QValue > 0.01.
        /// </summary>
        [Test]
        public static void WriteHighQValuePsmsFalse_OnlyLowQValueRowsInOutput()
        {
            string dir = MakeTempDir(nameof(WriteHighQValuePsmsFalse_OnlyLowQValueRowsInOutput));
            try
            {
                var protein = new CircularProtein("AFYRHTQESG", "CYC0003R00");
                string canonical = protein.BaseSequence;

                string mgfPath = WriteMgfForSequence(
                    Path.Combine(dir, "synthetic.mgf"), canonical,
                    maxMissedCleavages: 2, out _);

                string fastaPath = WriteFasta(
                    Path.Combine(dir, "synthetic.fasta"),
                    new[] { ("CYC0003R00", "CYC0003R00_SYNTH", canonical) });

                var searchParams = new CircularSearchParameters { WriteHighQValuePsms = false };

                RunTask(dir,
                    new List<DbForTask> { new DbForTask(fastaPath, isContaminant: false) },
                    mgfPath,
                    MakeTask(searchParams));

                var (headers, rows) = ReadOutputTsv(dir);

                foreach (var row in rows)
                {
                    string value = GetColumn(headers, row, "QValue");
                    if (double.TryParse(value,
                            System.Globalization.NumberStyles.Float,
                            System.Globalization.CultureInfo.InvariantCulture,
                            out double qValue))
                    {
                        Assert.That(qValue, Is.LessThanOrEqualTo(0.01),
                            $"Row with QValue={qValue} should not appear when " +
                            $"WriteHighQValuePsms=false.");
                    }
                }
            }
            finally { Directory.Delete(dir, recursive: true); }
        }
    }
}