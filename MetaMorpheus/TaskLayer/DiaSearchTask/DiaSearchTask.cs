// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.DiaSearch;
using MassSpectrometry;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskLayer
{
    /// <summary>
    /// MetaMorpheus task for DIA (Data-Independent Acquisition) spectral library search.
    /// 
    /// Orchestrates: load raw file → load spectral library → DiaEngine search → TSV output.
    /// Follows the same pattern as <see cref="SearchTask"/> for file loading, status reporting,
    /// and result writing.
    /// 
    /// Usage via CLI:
    ///   CMD.exe -d library.msp -s data.mzML -t DiaSearchTask.toml
    /// </summary>
    public class DiaSearchTask : MetaMorpheusTask
    {
        public MetaMorpheusDiaSearchParameters DiaSearchParameters { get; set; }

        public DiaSearchTask() : base(MyTask.Search) // reuse Search task type for now
        {
            CommonParameters = new CommonParameters();
            DiaSearchParameters = new MetaMorpheusDiaSearchParameters();
        }

        protected override MyTaskResults RunSpecific(
            string OutputFolder,
            List<DbForTask> dbFilenameList,
            List<string> currentRawFileList,
            string taskId,
            FileSpecificParameters[] fileSettingsList)
        {
            MyTaskResults = new MyTaskResults(this);
            var overallStopwatch = Stopwatch.StartNew();

            // ── Load spectral library ──────────────────────────────────────
            Status("Loading spectral library...", new List<string> { taskId });

            // Use the inherited LoadSpectralLibraries method (handles .msp files in dbFilenameList)
            // The user provides .msp files via the -d flag, same as SearchTask
            var spectralLibrary = LoadSpectralLibraries(taskId, dbFilenameList);

            if (spectralLibrary == null || spectralLibrary.Results == null || spectralLibrary.Results.Count == 0)
            {
                Warn("No spectral library found. DIA search requires a spectral library (.msp file). " +
                     "Include one in the database file list via the -d flag (e.g., -d library.msp).");
                return MyTaskResults;
            }

            List<LibrarySpectrum> librarySpectra = spectralLibrary.Results;
            Status($"Loaded {librarySpectra.Count:N0} library spectra",
                new List<string> { taskId });

            // ── Convert mzLib parameters ───────────────────────────────────
            var mzLibParams = DiaSearchParameters.ToMzLibParameters();

            // ── Prepare file manager ───────────────────────────────────────
            MyFileManager myFileManager = new MyFileManager(true); // dispose files after use

            // Start background loading of first file
            string firstFile = currentRawFileList[0];
            CommonParameters firstFileParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[0]);
            Task<MsDataFile> nextFileLoadingTask = new Task<MsDataFile>(
                () => myFileManager.LoadFile(firstFile, firstFileParams));
            nextFileLoadingTask.Start();

            // ── Aggregate results across all files ─────────────────────────
            var allResults = new List<(string FileName, List<MassSpectrometry.Dia.DiaSearchResult> Results)>();
            var diagnosticLines = new List<string>();

            for (int fileIndex = 0; fileIndex < currentRawFileList.Count; fileIndex++)
            {
                if (GlobalVariables.StopLoops) break;

                string rawFilePath = currentRawFileList[fileIndex];
                string rawFileName = Path.GetFileNameWithoutExtension(rawFilePath);
                var thisId = new List<string> { taskId, "Individual Spectra Files", rawFilePath };

                StartingDataFile(rawFilePath, thisId);
                Status("Loading spectra file...", thisId);

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(
                    CommonParameters, fileSettingsList[fileIndex]);

                // Wait for current file to finish loading
                nextFileLoadingTask.Wait();
                MsDataFile myMsDataFile = nextFileLoadingTask.Result;

                // Start background loading of next file
                if (fileIndex < currentRawFileList.Count - 1)
                {
                    int nextIndex = fileIndex + 1;
                    nextFileLoadingTask = new Task<MsDataFile>(
                        () => myFileManager.LoadFile(
                            currentRawFileList[nextIndex],
                            SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[nextIndex])));
                    nextFileLoadingTask.Start();
                }

                // ── Run DIA Engine ─────────────────────────────────────────
                Status("Running DIA search...", thisId);
                var fileStopwatch = Stopwatch.StartNew();

                MsDataScan[] scans = myMsDataFile.GetAllScansList().ToArray();

                using var diaEngine = new DiaEngine(
                    scans,
                    librarySpectra,
                    mzLibParams,
                    combinedParams,
                    this.FileSpecificParameters,
                    thisId);

                var engineResults = (DiaEngineResults)diaEngine.Run();
                fileStopwatch.Stop();

                myFileManager.DoneWithFile(rawFilePath);

                // ── Collect results ────────────────────────────────────────
                allResults.Add((rawFileName, engineResults.DiaResults));

                int targetCount = engineResults.TargetCount;
                int decoyCount = engineResults.DecoyCount;

                Status($"DIA search complete: {targetCount:N0} targets, {decoyCount:N0} decoys " +
                       $"in {fileStopwatch.Elapsed.TotalSeconds:F1}s", thisId);

                // ── Diagnostics ────────────────────────────────────────────
                if (DiaSearchParameters.WriteDiagnostics)
                {
                    diagnosticLines.Add(
                        $"{rawFileName}\t{scans.Length}\t{librarySpectra.Count}\t" +
                        $"{targetCount}\t{decoyCount}\t{fileStopwatch.Elapsed.TotalSeconds:F2}");
                }

                // ── Per-file TSV ───────────────────────────────────────────
                string perFileFolder = Path.Combine(OutputFolder, "Individual File Results");
                if (!Directory.Exists(perFileFolder))
                    Directory.CreateDirectory(perFileFolder);

                string perFileTsvPath = Path.Combine(perFileFolder, rawFileName + "_DiaResults.tsv");
                WriteDiaResultsTsv(perFileTsvPath, rawFileName, engineResults.DiaResults,
                    DiaSearchParameters.WriteDecoyResults);

                FinishedDataFile(rawFilePath, thisId);
                ReportProgress(new ProgressEventArgs(
                    (int)((fileIndex + 1.0) / currentRawFileList.Count * 100),
                    "DIA search complete!", thisId));
            }

            // ── Write aggregate TSV ────────────────────────────────────────
            if (allResults.Count > 0)
            {
                string aggregatePath = Path.Combine(OutputFolder, "AllDiaResults.tsv");
                WriteAggregateTsv(aggregatePath, allResults, DiaSearchParameters.WriteDecoyResults);
                ProseCreatedWhileRunning.AppendLine($"Aggregate DIA results written to {aggregatePath}");
            }

            // ── Write diagnostics ──────────────────────────────────────────
            if (DiaSearchParameters.WriteDiagnostics && diagnosticLines.Count > 0)
            {
                string diagPath = Path.Combine(OutputFolder, "DiaDiagnostics.tsv");
                WriteDiagnostics(diagPath, diagnosticLines, overallStopwatch.Elapsed);
            }

            overallStopwatch.Stop();

            int totalTargets = allResults.Sum(r => r.Results.Count(x => !x.IsDecoy));
            int totalDecoys = allResults.Sum(r => r.Results.Count(x => x.IsDecoy));
            ProseCreatedWhileRunning.AppendLine(
                $"DIA search complete: {totalTargets:N0} target identifications" +
                $" across {currentRawFileList.Count} file(s)" +
                $" in {overallStopwatch.Elapsed.TotalSeconds:F1}s");

            return MyTaskResults;
        }

        #region TSV Output

        /// <summary>Header row for DIA results TSV files.</summary>
        private static readonly string TsvHeader = string.Join("\t", new[]
        {
            "File Name",
            "Sequence",
            "Precursor m/z",
            "Charge",
            "Window ID",
            "Is Decoy",
            "Dot Product Score",
            "Spectral Angle Score",
            "Fragments Detected",
            "Fragments Queried",
            "Fragment Detection Rate",
            "Library RT",
            "RT Window Start",
            "RT Window End",
            "XIC Point Counts",
            "Extracted Intensities"
        });

        /// <summary>
        /// Writes DIA search results to a tab-separated file.
        /// </summary>
        private static void WriteDiaResultsTsv(
            string filePath,
            string rawFileName,
            List<MassSpectrometry.Dia.DiaSearchResult> results,
            bool writeDecoys)
        {
            using var writer = new StreamWriter(filePath, false, Encoding.UTF8);
            writer.WriteLine(TsvHeader);

            foreach (var result in results)
            {
                if (!writeDecoys && result.IsDecoy)
                    continue;

                writer.WriteLine(FormatResultRow(rawFileName, result));
            }
        }

        /// <summary>
        /// Writes aggregate results from multiple files into one TSV.
        /// </summary>
        private static void WriteAggregateTsv(
            string filePath,
            List<(string FileName, List<MassSpectrometry.Dia.DiaSearchResult> Results)> allResults,
            bool writeDecoys)
        {
            using var writer = new StreamWriter(filePath, false, Encoding.UTF8);
            writer.WriteLine(TsvHeader);

            foreach (var (fileName, results) in allResults)
            {
                foreach (var result in results)
                {
                    if (!writeDecoys && result.IsDecoy)
                        continue;

                    writer.WriteLine(FormatResultRow(fileName, result));
                }
            }
        }

        /// <summary>
        /// Formats a single DiaSearchResult as a TSV row.
        /// </summary>
        private static string FormatResultRow(string fileName, MassSpectrometry.Dia.DiaSearchResult r)
        {
            string xicCounts = r.XicPointCounts != null
                ? string.Join(";", r.XicPointCounts)
                : "";

            string extractedIntensities = r.ExtractedIntensities != null
                ? string.Join(";", r.ExtractedIntensities.Select(x => x.ToString("G6", CultureInfo.InvariantCulture)))
                : "";

            return string.Join("\t",
                fileName,
                r.Sequence,
                r.PrecursorMz.ToString("F4", CultureInfo.InvariantCulture),
                r.ChargeState.ToString(),
                r.WindowId.ToString(),
                r.IsDecoy ? "TRUE" : "FALSE",
                float.IsNaN(r.DotProductScore) ? "NaN" : r.DotProductScore.ToString("F4", CultureInfo.InvariantCulture),
                float.IsNaN(r.SpectralAngleScore) ? "NaN" : r.SpectralAngleScore.ToString("F4", CultureInfo.InvariantCulture),
                r.FragmentsDetected.ToString(),
                r.FragmentsQueried.ToString(),
                r.FragmentDetectionRate.ToString("F4", CultureInfo.InvariantCulture),
                r.LibraryRetentionTime.HasValue
                    ? r.LibraryRetentionTime.Value.ToString("F2", CultureInfo.InvariantCulture)
                    : "",
                r.RtWindowStart.ToString("F2", CultureInfo.InvariantCulture),
                r.RtWindowEnd.ToString("F2", CultureInfo.InvariantCulture),
                xicCounts,
                extractedIntensities);
        }

        /// <summary>
        /// Writes timing and throughput diagnostics.
        /// </summary>
        private static void WriteDiagnostics(
            string filePath,
            List<string> diagnosticLines,
            TimeSpan totalElapsed)
        {
            using var writer = new StreamWriter(filePath, false, Encoding.UTF8);
            writer.WriteLine("File\tScans\tLibrary Spectra\tTargets\tDecoys\tTime (s)");
            foreach (var line in diagnosticLines)
                writer.WriteLine(line);
            writer.WriteLine();
            writer.WriteLine($"Total elapsed: {totalElapsed.TotalSeconds:F2}s");
        }

        #endregion
    }
}