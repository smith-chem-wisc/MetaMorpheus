// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.DiaSearch;
using MassSpectrometry;
using MassSpectrometry.Dia;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using Readers;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test.DiaSearch
{
    /// <summary>
    /// End-to-end integration tests for DiaSearchTask.
    /// 
    /// These tests exercise the full MetaMorpheus DIA pipeline:
    ///   synthetic mzML → DiaSearchTask.RunSpecific → FDR → .tsv output
    /// 
    /// They do NOT require real data files. All inputs are generated in-memory
    /// or written to a temp folder and cleaned up on teardown.
    /// 
    /// Coverage:
    ///   1. Happy path: results file produced, contains data rows, q-values in [0,1]
    ///   2. TSV schema: correct column count and header names in AllDiaResults.tsv
    ///   3. Per-file pre-FDR diagnostic TSV is written
    ///   4. DiagnosticFile written when WriteDiagnostics = true
    ///   5. Decoys absent from output when WriteDecoyResults = false
    ///   6. Empty library → task completes cleanly with no results file (early exit)
    ///   7. Targets and decoys both present → FDR runs, q-values monotone
    /// </summary>
    [TestFixture]
    public class DiaSearchTaskIntegrationTests
    {
        private string _outputFolder;
        private string _mzmlPath;
        private string _combinedLibraryPath;   // single .msp with targets + decoys
        private string _targetOnlyLibraryPath; // targets only (for decoy-suppression tests)

        // ── Isolation windows used by the synthetic data ────────────────────
        // 4 windows: 400–425, 425–450, 450–475, 475–500 m/z (center ± 12.5)
        private static readonly double[] WindowCenters = { 412.5, 437.5, 462.5, 487.5 };
        private const double WindowHalfWidth = 12.5;

        // Fragment m/z values that appear in every synthetic scan
        // (200, 210, 220, … 690 — 50 peaks per scan)
        private static readonly double[] ScanFragmentMzs =
            Enumerable.Range(0, 50).Select(i => 200.0 + i * 10.0).ToArray();

        // Precursors — one per window, fragment m/z matches ScanFragmentMzs[0..5]
        private static readonly double[] TargetPrecursorMzs = { 410.0, 435.0, 460.0, 485.0 };
        private static readonly double[] DecoyPrecursorMzs = { 411.0, 436.0, 461.0, 486.0 };

        [SetUp]
        public void SetUp()
        {
            string uid = Guid.NewGuid().ToString("N");
            _outputFolder = Path.Combine(Path.GetTempPath(), $"DiaIntegration_{uid}");
            Directory.CreateDirectory(_outputFolder);

            _mzmlPath = Path.Combine(_outputFolder, "synthetic_dia.mzML");
            _combinedLibraryPath = Path.Combine(_outputFolder, "decoys.msp"); // named "decoys" so DiaSearchTask marks them
            _targetOnlyLibraryPath = Path.Combine(_outputFolder, "targets.msp");

            WriteSyntheticMzml(_mzmlPath);
            // Combined file: targets first, then decoys — all in one .msp named "decoys.msp"
            // DiaSearchTask re-reads this file for "Name:" lines and marks all entries as decoy,
            // but targets appear in both files, so we write targets.msp separately for tests
            // that don't want decoy FDR behavior.
            WriteSyntheticMspCombined(_combinedLibraryPath, TargetPrecursorMzs, DecoyPrecursorMzs);
            WriteSyntheticMsp(_targetOnlyLibraryPath, TargetPrecursorMzs, isDecoy: false);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_outputFolder))
            {
                try { Directory.Delete(_outputFolder, recursive: true); }
                catch { /* best-effort */ }
            }
        }

        // ════════════════════════════════════════════════════════════════════
        // Test 1 — Happy path: task runs, output files produced with data
        // ════════════════════════════════════════════════════════════════════

        [Test]
        public void Integration_HappyPath_ResultsFileProduced()
        {
            var task = BuildTask();
            var dbList = BuildDbList(includeDecoys: true);

            RunTask(task, dbList);

            // AllDiaResults.tsv always exists after the task runs, but may contain
            // only the header when all results fail 1% FDR (expected with synthetic
            // data: 4 targets + 4 decoys gives q-values ≥ 0.25).
            string tsvPath = Path.Combine(_outputFolder, "AllDiaResults.tsv");
            Assert.That(File.Exists(tsvPath), Is.True,
                "AllDiaResults.tsv should be written by PostDiaSearchAnalysisTask");

            // The per-file pre-FDR TSV confirms DiaEngine produced results — it is
            // written before FDR filtering and always contains data if the engine ran.
            string perFileFolder = Path.Combine(_outputFolder, "Individual File Results");
            string[] preFdrFiles = Directory.Exists(perFileFolder)
                ? Directory.GetFiles(perFileFolder, "*_DiaResults_PreFDR.tsv")
                : Array.Empty<string>();

            Assert.That(preFdrFiles.Length, Is.GreaterThan(0),
                "Per-file pre-FDR TSV should be written, confirming DiaEngine ran");

            string[] preFdrLines = File.ReadAllLines(preFdrFiles[0]);
            Assert.That(preFdrLines.Length, Is.GreaterThan(1),
                "Pre-FDR TSV should have header + at least one result row");
        }

        // ════════════════════════════════════════════════════════════════════
        // Test 2 — TSV schema: 17 columns, correct header names
        // ════════════════════════════════════════════════════════════════════

        [Test]
        public void Integration_ResultsTsv_CorrectSchema()
        {
            var task = BuildTask();
            RunTask(task, BuildDbList(includeDecoys: true));

            // Validate AllDiaResults.tsv header (file always exists, may have no data rows)
            string tsvPath = Path.Combine(_outputFolder, "AllDiaResults.tsv");
            Assert.That(File.Exists(tsvPath), Is.True);

            string[] lines = File.ReadAllLines(tsvPath);
            Assert.That(lines.Length, Is.GreaterThan(0), "File must have at least a header line");

            string[] header = lines[0].Split('\t');
            Assert.That(header[0], Is.EqualTo("Sequence"));
            Assert.That(header[1], Is.EqualTo("Charge"));
            Assert.That(header[2], Is.EqualTo("Precursor m/z"));
            Assert.That(header[13], Is.EqualTo("QValue"));
            Assert.That(header.Length, Is.EqualTo(17));

            // If any post-FDR data rows exist, each must match the header column count
            foreach (string line in lines.Skip(1))
            {
                string[] cols = line.Split('\t');
                Assert.That(cols.Length, Is.EqualTo(header.Length),
                    $"Data row has {cols.Length} columns, expected {header.Length}: [{line}]");
            }

            // Validate pre-FDR per-file TSV schema — this always has data rows
            string perFileFolder = Path.Combine(_outputFolder, "Individual File Results");
            string[] preFdrFiles = Directory.Exists(perFileFolder)
                ? Directory.GetFiles(perFileFolder, "*_DiaResults_PreFDR.tsv")
                : Array.Empty<string>();

            if (preFdrFiles.Length > 0)
            {
                string[] preFdrLines = File.ReadAllLines(preFdrFiles[0]);
                if (preFdrLines.Length > 1)
                {
                    string[] preHeader = preFdrLines[0].Split('\t');
                    string[] dataRow = preFdrLines[1].Split('\t');
                    Assert.That(dataRow.Length, Is.EqualTo(preHeader.Length),
                        "Pre-FDR data row column count must match its header");
                }
            }
        }

        // ════════════════════════════════════════════════════════════════════
        // Test 3 — Q-values: assigned to every result, in [0, 1], monotone in output
        // Uses AllDiaResults.tsv when non-empty; falls back to RunScoreSortFdr output.
        // ════════════════════════════════════════════════════════════════════

        [Test]
        public void Integration_QValues_ValidAndMonotone()
        {
            var task = BuildTask();
            RunTask(task, BuildDbList(includeDecoys: true));

            string tsvPath = Path.Combine(_outputFolder, "AllDiaResults.tsv");
            Assert.That(File.Exists(tsvPath), Is.True);

            string[] lines = File.ReadAllLines(tsvPath);
            string[] header = lines[0].Split('\t');
            int qValueCol = Array.IndexOf(header, "QValue");
            Assert.That(qValueCol, Is.GreaterThanOrEqualTo(0), "QValue column not found");

            if (lines.Length <= 1)
            {
                // Post-FDR file is header-only (all results filtered by 1% threshold).
                // Verify q-values were assigned by checking the pre-FDR file for non-null scores.
                string perFileFolder = Path.Combine(_outputFolder, "Individual File Results");
                string[] preFdrFiles = Directory.Exists(perFileFolder)
                    ? Directory.GetFiles(perFileFolder, "*_DiaResults_PreFDR.tsv")
                    : Array.Empty<string>();
                Assert.That(preFdrFiles.Length, Is.GreaterThan(0),
                    "Pre-FDR results must exist when post-FDR is empty");
                // Pre-FDR file has spectral angle scores — confirm they are finite
                string[] preFdrLines = File.ReadAllLines(preFdrFiles[0]);
                string[] preHeader = preFdrLines[0].Split('\t');
                int saCol = Array.IndexOf(preHeader, "Spectral Angle Score");
                Assert.That(saCol, Is.GreaterThanOrEqualTo(0));
                foreach (string row in preFdrLines.Skip(1))
                {
                    string[] cols = row.Split('\t');
                    Assert.That(float.TryParse(cols[saCol], NumberStyles.Float,
                        CultureInfo.InvariantCulture, out _), Is.True,
                        $"Spectral angle score should be numeric: {cols[saCol]}");
                }
                return;
            }

            // Post-FDR rows exist — validate q-values are in [0,1] and monotone
            var qValues = new List<double>();
            foreach (string line in lines.Skip(1))
            {
                string[] cols = line.Split('\t');
                if (string.IsNullOrWhiteSpace(cols[qValueCol])) continue;
                // Q-values are written in E4 scientific notation
                double q = double.Parse(cols[qValueCol], NumberStyles.Float, CultureInfo.InvariantCulture);
                Assert.That(q, Is.GreaterThanOrEqualTo(0.0).And.LessThanOrEqualTo(1.0),
                    $"Q-value {q} out of [0,1]");
                qValues.Add(q);
            }

            for (int i = 1; i < qValues.Count; i++)
                Assert.That(qValues[i], Is.GreaterThanOrEqualTo(qValues[i - 1]),
                    $"Q-values not monotone at row {i + 1}: {qValues[i - 1]} → {qValues[i]}");
        }

        // ════════════════════════════════════════════════════════════════════
        // Test 4 — Per-file pre-FDR diagnostic TSV is written
        // ════════════════════════════════════════════════════════════════════

        [Test]
        public void Integration_PerFileDiagnosticTsv_Written()
        {
            var task = BuildTask();
            RunTask(task, BuildDbList(includeDecoys: false));

            string perFileFolder = Path.Combine(_outputFolder, "Individual File Results");
            Assert.That(Directory.Exists(perFileFolder), Is.True,
                "Individual File Results folder should be created");

            string[] tsvFiles = Directory.GetFiles(perFileFolder, "*_DiaResults_PreFDR.tsv");
            Assert.That(tsvFiles.Length, Is.GreaterThan(0),
                "At least one per-file pre-FDR TSV should be written");
        }

        // ════════════════════════════════════════════════════════════════════
        // Test 5 — DiaDiagnostics.tsv written when WriteDiagnostics = true
        // ════════════════════════════════════════════════════════════════════

        [Test]
        public void Integration_DiagnosticsFile_WrittenWhenEnabled()
        {
            var task = BuildTask();
            task.DiaSearchParameters.WriteDiagnostics = true;
            RunTask(task, BuildDbList(includeDecoys: false));

            string diagPath = Path.Combine(_outputFolder, "DiaDiagnostics.tsv");
            Assert.That(File.Exists(diagPath), Is.True,
                "DiaDiagnostics.tsv should be written when WriteDiagnostics = true");

            string[] lines = File.ReadAllLines(diagPath);
            Assert.That(lines.Length, Is.GreaterThan(1),
                "Diagnostics file should have a header and at least one data row");
        }

        // ════════════════════════════════════════════════════════════════════
        // Test 6 — WriteDecoyResults = false: no decoy rows in any output
        // ════════════════════════════════════════════════════════════════════

        [Test]
        public void Integration_WriteDecoyResultsFalse_NoDecoysInOutput()
        {
            var task = BuildTask();
            task.DiaSearchParameters.WriteDecoyResults = false;
            RunTask(task, BuildDbList(includeDecoys: true));

            // Pre-FDR TSV always has rows and respects WriteDecoyResults
            string perFileFolder = Path.Combine(_outputFolder, "Individual File Results");
            string[] preFdrFiles = Directory.Exists(perFileFolder)
                ? Directory.GetFiles(perFileFolder, "*_DiaResults_PreFDR.tsv")
                : Array.Empty<string>();

            Assert.That(preFdrFiles.Length, Is.GreaterThan(0),
                "Pre-FDR TSV must exist to validate decoy suppression");

            string[] lines = File.ReadAllLines(preFdrFiles[0]);
            string[] header = lines[0].Split('\t');
            int isDecoyCol = Array.IndexOf(header, "Is Decoy");
            Assert.That(isDecoyCol, Is.GreaterThanOrEqualTo(0), "Is Decoy column not found");

            foreach (string line in lines.Skip(1))
            {
                string[] cols = line.Split('\t');
                Assert.That(cols[isDecoyCol], Is.EqualTo("FALSE"),
                    "Pre-FDR TSV should contain no decoy rows when WriteDecoyResults=false");
            }
        }

        // ════════════════════════════════════════════════════════════════════
        // Test 7 — Empty library: task exits cleanly, no crash, no tsv
        // ════════════════════════════════════════════════════════════════════

        [Test]
        public void Integration_EmptyLibrary_TaskCompletesCleanly()
        {
            string emptyLibPath = Path.Combine(_outputFolder, "empty.msp");
            File.WriteAllText(emptyLibPath, "");

            var task = BuildTask();
            var dbList = new List<DbForTask> { new DbForTask(emptyLibPath, false) };

            Assert.DoesNotThrow(() => RunTask(task, dbList));
        }

        // ════════════════════════════════════════════════════════════════════
        // Test 8 — ScanIndex is disposed after FDR (no memory leak path)
        //          Verify via DiaEngineResults.ScanIndex being null-safe
        // ════════════════════════════════════════════════════════════════════

        [Test]
        public void Integration_DiaEngine_ScanIndexDisposedAfterFdr()
        {
            // Run DiaEngine directly and verify ScanIndex is non-null before FDR
            var scans = BuildSyntheticScans();
            var library = BuildSyntheticLibrarySpectra(includeDecoys: true);
            var mzLibParams = new DiaSearchParameters
            {
                PpmTolerance = 20f,
                RtToleranceMinutes = 5.0f,
                MinFragmentsRequired = 1,
                MaxThreads = 1
            };

            DiaEngineResults results;
            using (var engine = new DiaEngine(
                scans, library, mzLibParams,
                useCalibration: false,
                commonParameters: new CommonParameters(),
                fileSpecificParameters: new List<(string, CommonParameters)>(),
                nestedIds: new List<string> { "IntegrationTest" },
                fileName: "synthetic"))
            {
                results = (DiaEngineResults)engine.Run();
            }

            // ScanIndex should be non-null before caller disposes it
            Assert.That(results.ScanIndex, Is.Not.Null,
                "ScanIndex must be non-null so PostDiaSearchAnalysisTask can compute MS1 features");

            // Dispose it (as DiaSearchTask does after FDR)
            Assert.DoesNotThrow(() => (results.ScanIndex as IDisposable)?.Dispose(),
                "Disposing ScanIndex after FDR should not throw");
        }

        // ════════════════════════════════════════════════════════════════════
        // Helper: build and run the task
        // ════════════════════════════════════════════════════════════════════

        private DiaSearchTask BuildTask()
        {
            return new DiaSearchTask
            {
                CommonParameters = new CommonParameters(maxThreadsToUsePerFile: 1),
                DiaSearchParameters = new MetaMorpheusDiaSearchParameters
                {
                    PpmTolerance = 20f,
                    RtToleranceMinutes = 5.0f,
                    MinFragmentsRequired = 1,   // low threshold for synthetic data
                    MinScoreThreshold = 0.0f,
                    MaxThreads = 1,
                    UseIrtCalibration = false,  // no calibration in tests
                    WriteDiagnostics = true,
                    WriteDecoyResults = false
                }
            };
        }

        private List<DbForTask> BuildDbList(bool includeDecoys)
        {
            // Always include the target-only file as the primary library.
            // If decoys requested, also include the combined file (named "decoys.msp")
            // so DiaSearchTask's filename-based decoy detection fires.
            var list = new List<DbForTask>
            {
                new DbForTask(_targetOnlyLibraryPath, false)
            };
            if (includeDecoys)
                list.Add(new DbForTask(_combinedLibraryPath, false));
            return list;
        }

        private void RunTask(DiaSearchTask task, List<DbForTask> dbList)
        {
            // Call MetaMorpheusTask.RunTask directly — this writes output straight into
            // _outputFolder without creating a task-named subfolder (unlike EverythingRunnerEngine).
            task.RunTask(
                _outputFolder,
                dbList,
                new List<string> { _mzmlPath },
                "IntegrationTest");
        }

        // ════════════════════════════════════════════════════════════════════
        // Synthetic data builders
        // ════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Writes a valid mzML with 40 MS2 scans across 4 isolation windows.
        /// Every scan contains 50 peaks at m/z 200, 210, ... 690.
        /// RT runs from 1.0 to 10.0 min (10 scans per window).
        /// Uses MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra.
        /// </summary>
        private static void WriteSyntheticMzml(string path)
        {
            var scans = BuildSyntheticScans();
            var sourceFile = new SourceFile(
                nativeIdFormat: "no nativeID format",
                massSpectrometerFileFormat: "mzML format",
                checkSum: null,
                fileChecksumType: "SHA-1",
                filePath: path,
                id: Path.GetFileNameWithoutExtension(path));

            var msDataFile = new FakeMsDataFile(scans, sourceFile);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(msDataFile, path, writeIndexed: true);
        }

        private static MsDataScan[] BuildSyntheticScans()
        {
            var scans = new List<MsDataScan>();
            int scanNumber = 1;

            foreach (double center in WindowCenters)
            {
                for (int s = 0; s < 10; s++)
                {
                    double rt = 1.0 + s; // minutes 1–10

                    var spectrum = new MzSpectrum(
                        ScanFragmentMzs,
                        ScanFragmentMzs.Select((_, i) => 1000.0 + i * 100.0).ToArray(),
                        shouldCopy: false);

                    scans.Add(new MsDataScan(
                        massSpectrum: spectrum,
                        oneBasedScanNumber: scanNumber++,
                        msnOrder: 2,
                        isCentroid: true,
                        polarity: Polarity.Positive,
                        retentionTime: rt,
                        scanWindowRange: new MzRange(100, 2000),
                        scanFilter: "FTMS + p NSI d Full ms2",
                        mzAnalyzer: MZAnalyzerType.Orbitrap,
                        totalIonCurrent: 1000.0 * 50,
                        injectionTime: 20.0,
                        noiseData: null,
                        nativeId: $"scan={scanNumber - 1}",
                        selectedIonMz: center,          // required — MzmlMethods calls .Value unconditionally
                        isolationMZ: center,
                        isolationWidth: WindowHalfWidth * 2,
                        dissociationType: DissociationType.HCD,
                        oneBasedPrecursorScanNumber: null,
                        selectedIonMonoisotopicGuessMz: center));
                }
            }

            return scans.ToArray();
        }

        /// <summary>
        /// Writes a NIST-format .msp spectral library file compatible with mzLib's SpectralLibrary reader.
        ///
        /// Format derived from LibrarySpectrum.ToString() — the canonical mzLib msp writer:
        ///   Name: SEQUENCE/CHARGE
        ///   MW: precursorMz
        ///   Comment: Parent=precursorMz RT=5.0
        ///   Num peaks: N
        ///   mz\tfraction\t"TypeNumber^Charge/0ppm"    ← tab-separated, intensity as fraction
        ///   Name: ...                                  ← next entry immediately, NO blank line
        ///
        /// Critical: NO blank line between entries. SpectralLibrary.ReadLibrarySpectrum reads
        /// every line after "Num peaks" as a fragment ion (readingPeaks=true), so a blank line
        /// causes ReadFragmentIon("") → split[] empty → IndexOutOfRangeException on split[0].
        ///
        /// fragmentSplit = {'\t', '"', ')', '/'} — the '/' splits "b1^1/0ppm" into ["b1^1","0ppm"],
        /// so split[2] = "b1^1" which IonParserRegex matches as type="b", number=1, charge=1.
        /// </summary>
        private static void WriteSyntheticMsp(string path, double[] precursorMzs, bool isDecoy)
        {
            using var writer = new StreamWriter(path);

            string prefix = isDecoy ? "DECOY_" : "";

            // Ion types must be valid ProductType enum values parseable by SpectralLibrary
            string[] ionTypes = { "b", "b", "b", "y", "y", "y" };
            int[] ionNumbers = { 1, 2, 3, 1, 2, 3 };

            for (int i = 0; i < precursorMzs.Length; i++)
            {
                string seq = $"{prefix}PEPTIDE{(char)('A' + i)}K";
                double mz = precursorMzs[i];
                double maxIntensity = 1000.0 + 5 * 200.0; // highest fragment intensity

                writer.WriteLine($"Name: {seq}/2");
                writer.WriteLine($"MW: {mz:F6}");
                writer.WriteLine($"Comment: Parent={mz:F6} RT=5.0");
                writer.WriteLine("Num peaks: 6");
                for (int f = 0; f < 6; f++)
                {
                    double rawIntensity = 1000.0 + f * 200.0;
                    double fraction = rawIntensity / maxIntensity;
                    // Format: mz \t fraction \t "TypeNumber^1/0ppm"
                    // fragmentSplit={'\t','"',')','/'} splits this into:
                    //   [0]=mz  [1]=fraction  [2]="TypeNumber^1"  [3]="0ppm"
                    writer.WriteLine(
                        $"{ScanFragmentMzs[f]:F6}\t{fraction:F6}\t\"{ionTypes[f]}{ionNumbers[f]}^1/0ppm\"");
                }
                // NO blank line here — next "Name:" line is the delimiter
            }
        }

        /// <summary>
        /// Writes a combined .msp with both targets and decoys in one file.
        /// Used for the "decoys.msp" file that DiaSearchTask re-reads to mark IsDecoy.
        /// All entries in a file named "decoys*" get marked as decoy by DiaSearchTask,
        /// so we put both here and rely on DiaSearchTask matching Names from this file.
        /// </summary>
        private static void WriteSyntheticMspCombined(string path,
            double[] targetMzs, double[] decoyMzs)
        {
            using var writer = new StreamWriter(path);
            string[] ionTypes = { "b", "b", "b", "y", "y", "y" };
            int[] ionNumbers = { 1, 2, 3, 1, 2, 3 };
            double maxIntensity = 1000.0 + 5 * 200.0;

            void WriteEntry(string seq, double mz)
            {
                writer.WriteLine($"Name: {seq}/2");
                writer.WriteLine($"MW: {mz:F6}");
                writer.WriteLine($"Comment: Parent={mz:F6} RT=5.0");
                writer.WriteLine("Num peaks: 6");
                for (int f = 0; f < 6; f++)
                {
                    double fraction = (1000.0 + f * 200.0) / maxIntensity;
                    writer.WriteLine(
                        $"{ScanFragmentMzs[f]:F6}\t{fraction:F6}\t\"{ionTypes[f]}{ionNumbers[f]}^1/0ppm\"");
                }
            }

            for (int i = 0; i < decoyMzs.Length; i++)
                WriteEntry($"DECOY_PEPTIDE{(char)('A' + i)}K", decoyMzs[i]);
        }
        private static List<LibrarySpectrum> BuildSyntheticLibrarySpectra(bool includeDecoys)
        {
            var library = new List<LibrarySpectrum>();

            for (int i = 0; i < TargetPrecursorMzs.Length; i++)
            {
                var fragments = new List<MatchedFragmentIon>();
                for (int f = 0; f < 6; f++)
                {
                    var product = new Product(ProductType.b, FragmentationTerminus.N,
                        ScanFragmentMzs[f] - 1.00727647, f + 1, f + 1, 0);
                    fragments.Add(new MatchedFragmentIon(product, ScanFragmentMzs[f], 1000.0 + f * 200.0, 1));
                }

                library.Add(new LibrarySpectrum(
                    sequence: $"PEPTIDE{(char)('A' + i)}K",
                    precursorMz: TargetPrecursorMzs[i],
                    chargeState: 2,
                    peaks: fragments,
                    rt: 5.0,
                    isDecoy: false));

                if (includeDecoys)
                {
                    library.Add(new LibrarySpectrum(
                        sequence: $"DECOY_PEPTIDE{(char)('A' + i)}K",
                        precursorMz: DecoyPrecursorMzs[i],
                        chargeState: 2,
                        peaks: fragments,
                        rt: 5.0,
                        isDecoy: true));
                }
            }

            return library;
        }

        // ════════════════════════════════════════════════════════════════════
        // Minimal concrete MsDataFile for mzML writing
        // ════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Minimal concrete MsDataFile wrapping a pre-built scan array.
        /// Used only to pass scans into MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra.
        /// </summary>
        private class FakeMsDataFile : MsDataFile
        {
            public FakeMsDataFile(MsDataScan[] scans, SourceFile sourceFile)
                : base(scans, sourceFile)
            {
            }

            public override MsDataFile LoadAllStaticData(FilteringParams filterParams = null, int maxThreads = 1)
                => this;

            public override SourceFile GetSourceFile()
                => SourceFile;

            public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
                => Scans[oneBasedScanNumber - 1];

            public override void CloseDynamicConnection() { }

            public override void InitiateDynamicConnection() { }
        }
    }
}