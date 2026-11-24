using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.Truncation;
using MassSpectrometry;
using Nett;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;

namespace TaskLayer
{
    /// <summary>
    /// The truncation search a <see cref="SearchTask"/> runs on its own matches when <see cref="CommonParameters.AddTruncations"/>
    /// is set. It identifies N- and C-terminally truncated proteoforms by re-searching the MS2 scans against a fragment index
    /// built from the proteoforms the search found. See docs/Truncation-Search.md for the full three-pass design (index + dual
    /// single-series scoring, terminus-directed chopping, pooled FDR/PEP).
    ///
    /// <see cref="RunTruncationSearch"/> wires the engine pieces together: dedup the search's matches to proteoform-level
    /// parents, add reverse decoys, run Pass 2 dual single-series scoring per file, chop the winners (Pass 3), pool with the
    /// inherited intact matches, run pooled FDR/PEP, and write AllTruncatedPSMs/AllTruncatedProteoforms.psmtsv.
    /// </summary>
    public partial class SearchTask
    {
        public const string TruncatedPsmsFileName = "AllTruncatedPSMs.psmtsv";
        public const string TruncatedProteoformsFileName = "AllTruncatedProteoforms.psmtsv";

        /// <summary>
        /// Parent inclusion q-value threshold (decision #3). Deliberately permissive (10% FDR) to cast a wide net of parents;
        /// the truncation search's own FDR analysis polices the final output.
        /// </summary>
        public const double TruncationParentQValueThreshold = 0.10;

        /// <summary>
        /// Notch acceptor used when chopping a parent down to the precursor mass in Pass 3 (decisions #9, #15), which also sets
        /// the FDR notch count of the pooled analysis: the top-down ThreeMM notches, so precursor isotope errors do not block a chop.
        /// </summary>
        public const MassDiffAcceptorType TruncationChopAcceptorType = MassDiffAcceptorType.ThreeMM;

        /// <summary>
        /// Environment variable that turns on the truncation benchmarking instrumentation: set it to the path of a rolling,
        /// append-only <c>perf_log.tsv</c>, and each run appends one metrics+timing row and writes the CandidateRanks.tsv side file.
        /// </summary>
        public const string TruncationPerfLogPathEnvironmentVariable = "MetaMorpheusTruncationPerfLog";

        /// <summary>
        /// Path to the perf log (docs/Truncation-Search.md, Benchmarking hook), set directly by tests; when null the
        /// <see cref="TruncationPerfLogPathEnvironmentVariable"/> environment variable is read. Benchmarking
        /// instrumentation, not a search setting, so it is not in the TOML.
        /// </summary>
        [TomlIgnore]
        internal string TruncationPerfLogPath { get; set; }

        // Set during parent building; reported in the perf log (did PEP, not notch q, drive the #3 filter?).
        private bool _pepQWasUsedForParents;

        /// <summary>
        /// Runs the truncation search on this search's resolved, FDR'd matches and writes its files into this search's folder
        /// under its task id. The matches keep the matched Protein + DigestionParams, which Pass 3 chopping and decoy
        /// generation need.
        /// </summary>
        internal void RunTruncationSearch(List<SpectralMatch> searchMatches, string OutputFolder, List<string> currentRawFileList,
            string taskId, FileSpecificParameters[] fileSettingsList, List<Modification> fixedModifications)
        {
            _pepQWasUsedForParents = false; // the task object can be re-run (GUI), so do not carry the flag over
            var wallStopwatch = System.Diagnostics.Stopwatch.StartNew();
            string perfLogPath = TruncationPerfLogPath ?? Environment.GetEnvironmentVariable(TruncationPerfLogPathEnvironmentVariable);

            MassDiffAcceptor chopAcceptor = GetChopAcceptor(CommonParameters.PrecursorMassTolerance);

            // Per-file lookup of the search's matches, used to inherit intact matches as full-length forms (#4a).
            ILookup<string, SpectralMatch> pass1ByFile = searchMatches.Where(p => p != null).ToLookup(p => p.FullFilePath);

            List<TruncationParent> parents = BuildTruncationParents(searchMatches);

            if (parents.Count == 0)
            {
                Warn("The truncation search found no proteoforms to seed it; writing empty result files.");
                WriteTruncationOutputs(new List<SpectralMatch>(), OutputFolder, taskId);
                return;
            }

            // 2. Augment with reverse decoys (decision #14).
            int targetAndPass1DecoyParents = parents.Count;
            parents = TruncationParentBuilder.AddReverseDecoys(parents);
            int decoysGenerated = parents.Count - targetAndPass1DecoyParents;

            var pooled = new List<SpectralMatch>();
            // Diagnostic: winning parent's 0-based rank in each scan's index match-count ordering (validates
            // the engine's MaxCandidatesToScore cap). Written to CandidateRanks.tsv.
            var winnerRankRows = new List<string>();

            // Perf accumulators (docs/Truncation-Search.md, Benchmarking hook); parent/oversize counts are identical per file (shared list).
            var pass3Timings = new TruncationTimings();
            double pass2IndexSeconds = 0, pass2ScoringSeconds = 0;
            int totalMs2Scans = 0, indexedParents = 0, oversizeExcluded = 0;

            var myFileManager = new MyFileManager(true);
            for (int i = 0; i < currentRawFileList.Count; i++)
            {
                // Load one file's MS2 scans at a time, so only the file being searched is held in memory.
                string rawFilePath = currentRawFileList[i];
                CommonParameters fileParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[i], currentRawFileList[i]);

                MsDataFile dataFile = myFileManager.LoadFile(rawFilePath, fileParams);
                Ms2ScanWithSpecificMass[] scans = GetMs2Scans(dataFile, rawFilePath, fileParams)
                    .OrderBy(s => s.OneBasedScanNumber).ToArray();
                myFileManager.DoneWithFile(rawFilePath);

                // Every Pass 1 PSM on each scan, best first: a chimeric scan can carry several confident
                // intact matches, one per precursor, and the engine matches each precursor entry to its own (#4a).
                Dictionary<int, IReadOnlyList<SpectralMatch>> pass1PsmsByScan = pass1ByFile[rawFilePath]
                    .GroupBy(p => p.ScanNumber)
                    .ToDictionary(g => g.Key, g => (IReadOnlyList<SpectralMatch>)g.OrderByDescending(p => p).ToList());

                // Pass 3 chop acceptor at this file's precursor tolerance, the same one Pass 2 uses (#9, #15).
                MassDiffAcceptor fileChopAcceptor = GetChopAcceptor(fileParams.PrecursorMassTolerance);

                var file = (Raw: rawFilePath, FileParams: fileParams, Scans: scans);
                var engine = new TruncationSearchEngine(parents, file.Scans, file.FileParams,
                    new TruncationAcceptor(file.FileParams.PrecursorMassTolerance), pass1PsmsByScan,
                    TruncationSearchEngine.DefaultMaxFragmentSize);
                List<TruncationParentSelection> selections = engine.Run();
                foreach (string warning in engine.Warnings)
                {
                    Warn(warning);
                }

                // Distinct raw MS2 scans (scans[] is precursor-expanded: one entry per candidate precursor).
                totalMs2Scans += file.Scans.Select(s => s.OneBasedScanNumber).Distinct().Count();
                pass2IndexSeconds += engine.IndexBuildSeconds;
                pass2ScoringSeconds += engine.ScoringSeconds;
                indexedParents = engine.IndexedParentCount;
                oversizeExcluded = engine.ExcludedOversizedParentCount;

                var truncationPsms = new List<TruncationPsm>();
                foreach (TruncationParentSelection selection in selections)
                {
                    if (selection.Outcome == TruncationScanOutcome.Winner)
                    {
                        List<TruncationPsm> psms = TruncationPass3.ScoreTruncations(selection, file.FileParams,
                            fileChopAcceptor, pass3Timings, fixedModifications);
                        truncationPsms.AddRange(psms);
                        winnerRankRows.Add(string.Join("\t", Path.GetFileNameWithoutExtension(file.Raw),
                            selection.Scan.OneBasedScanNumber, selection.CandidateRank, selection.CandidatePoolSize,
                            selection.Score.ToString("F3", System.Globalization.CultureInfo.InvariantCulture),
                            psms.Count > 0 ? 1 : 0));
                    }
                    else if (selection.Outcome == TruncationScanOutcome.IntactInherited && selection.IntactMatch != null)
                    {
                        pooled.Add(TruncationPass3.InheritAsFullLength(selection.IntactMatch, selection.Scan, selection.ScanIndex, file.FileParams));
                    }
                }

                // One truncation PSM per scan, tied parents kept as ambiguity (#10), before pooling.
                pooled.AddRange(TruncationPass3.CollapseDuplicateTruncations(truncationPsms));
            }

            // 4. Pooled FDR + PEP over (full-length + truncation) PSMs (decisions #15, #18). An empty pool
            //    (no winners and no inherited intact matches) skips FDR and writes header-only files.
            var fdrStopwatch = System.Diagnostics.Stopwatch.StartNew();
            List<SpectralMatch> withFdr = pooled.Count == 0
                ? pooled
                : TruncationFdr.RunPooledFdr(pooled, CommonParameters,
                    chopAcceptor.NumNotches, FileSpecificParameters, taskId, OutputFolder);
            double fdrPepSeconds = fdrStopwatch.Elapsed.TotalSeconds;

            // 5. Write the two result files (decisions #16, #17). No q/PEP cutoff; decoys/contaminants as the search writes them.
            WriteTruncationOutputs(withFdr, OutputFolder, taskId);

            // Diagnostic side file: per-scan winning-parent rank in the index match-count ordering.
            // Opt-in (gated on the same perf-log path that enables benchmarking), registered via
            // FinishedWritingFile, and best-effort so a diagnostic-write failure can't abort an otherwise
            // successful search (#14).
            if (!string.IsNullOrWhiteSpace(perfLogPath))
            {
                try
                {
                    string ranksPath = Path.Combine(OutputFolder, "CandidateRanks.tsv");
                    File.WriteAllLines(ranksPath,
                        new[] { "File\tScanNumber\tCandidateRank\tCandidatePoolSize\tScore\tProducedTruncation" }.Concat(winnerRankRows));
                    FinishedWritingFile(ranksPath, new List<string> { taskId });
                }
                catch (Exception ex)
                {
                    Warn($"Could not write diagnostic CandidateRanks.tsv: {ex.Message}");
                }

                // 6. Perf-log row (docs/Truncation-Search.md, Benchmarking hook).
                AppendTruncationPerfLog(perfLogPath, withFdr, OutputFolder, currentRawFileList.Count, totalMs2Scans,
                    wallStopwatch.Elapsed.TotalSeconds, pass2IndexSeconds, pass2ScoringSeconds, pass3Timings,
                    fdrPepSeconds, indexedParents, oversizeExcluded, decoysGenerated);
            }
        }

        /// <summary>Computes the per-run metrics from the written set and appends one perf_log.tsv row.</summary>
        private void AppendTruncationPerfLog(string perfLogPath, List<SpectralMatch> withFdr, string outputFolder, int nRawFiles, int totalMs2Scans,
            double wallSeconds, double pass2IndexSeconds, double pass2ScoringSeconds, TruncationTimings pass3Timings,
            double fdrPepSeconds, int indexedParents, int oversizeExcluded, int decoysGenerated)
        {
            List<SpectralMatch> rows = withFdr.Where(p => p != null).ToList();
            int Truncations(string kind) => rows.Count(p => TruncationDescriptionOf(p).Contains(kind));
            int targetPsmsAtQ(double q) => rows.Count(p => !p.IsDecoy && (p.GetFdrInfo(false)?.QValue ?? 2) <= q);
            int targetProteoformsAtQ(double q) => rows.Where(p => !p.IsDecoy && (p.GetFdrInfo(true)?.QValueNotch ?? 2) <= q)
                .Select(p => p.FullSequence).Distinct().Count();

            var metrics = new TruncationPerfMetrics
            {
                OutputFolder = outputFolder,
                NRawFiles = nRawFiles,
                TotalMs2Scans = totalMs2Scans,
                TaskWallSeconds = wallSeconds,
                NPsmsEmitted = rows.Count,
                NProteoformsEmitted = rows.Select(p => p.FullSequence).Distinct().Count(),
                NPsmsQ01 = targetPsmsAtQ(0.01),
                NProteoformsQ01 = targetProteoformsAtQ(0.01),
                NPsmsQ05 = targetPsmsAtQ(0.05),
                NProteoformsQ05 = targetProteoformsAtQ(0.05),
                NTruncationsTotal = Truncations("truncation"),
                NTruncationsNterm = Truncations(TruncationPass3.NTerminalTruncation),
                NTruncationsCterm = Truncations(TruncationPass3.CTerminalTruncation),
                NIntactInherited = Truncations(TruncationPass3.FullLength),
                NParentsIndexed = indexedParents,
                NParentsOversizeExcluded = oversizeExcluded,
                PepQWasUsed = _pepQWasUsedForParents,
                NDecoysGenerated = decoysGenerated,
                Pass2IndexSeconds = pass2IndexSeconds,
                Pass2ScoringSeconds = pass2ScoringSeconds,
                Pass3ChoppingSeconds = pass3Timings.ChoppingSeconds,
                Pass3ScoringSeconds = pass3Timings.ScoringSeconds,
                FdrPepSeconds = fdrPepSeconds
            };

            // Run metadata from the parent run-folder name (<date>_<phase>_<datasetTag>_<runLabel>).
            string runFolder = Directory.GetParent(outputFolder)?.Name ?? "";
            (metrics.Phase, metrics.DatasetTag, metrics.RunLabel) = PerfLogger.ParseRunFolderName(runFolder);

            try
            {
                PerfLogger.Append(perfLogPath, metrics);
            }
            catch (Exception ex)
            {
                Warn("Failed to append perf-log row: " + ex.Message);
            }
        }

        private static string TruncationDescriptionOf(SpectralMatch psm)
        {
            var pwsm = psm.BestMatchingBioPolymersWithSetMods.FirstOrDefault()?.SpecificBioPolymer as PeptideWithSetModifications;
            return pwsm?.Description ?? string.Empty;
        }

        /// <summary>
        /// Writes AllTruncatedPSMs.psmtsv and AllTruncatedProteoforms.psmtsv (always, even when empty). No q/PEP cutoff; decoy
        /// and contaminant rows and the Essential Sequence mods follow this search's WriteDecoys / WriteContaminants /
        /// ModsToWriteSelection (#17), so the truncation files join to the search's AllPSMs/AllProteoforms.
        /// </summary>
        private void WriteTruncationOutputs(List<SpectralMatch> psms, string outputFolder, string taskId)
        {
            List<SpectralMatch> toWrite = TruncationRowsToWrite(psms, SearchParameters);

            string psmsPath = Path.Combine(outputFolder, TruncatedPsmsFileName);
            string proteoformsPath = Path.Combine(outputFolder, TruncatedProteoformsFileName);
            TruncationOutput.WritePsms(toWrite, psmsPath, SearchParameters.ModsToWriteSelection);
            TruncationOutput.WriteProteoforms(toWrite, proteoformsPath, SearchParameters.ModsToWriteSelection);
            FinishedWritingFile(psmsPath, new List<string> { taskId });
            FinishedWritingFile(proteoformsPath, new List<string> { taskId });
        }

        /// <summary>The PSMs that reach the truncation output files: decoys and contaminants only when WriteDecoys / WriteContaminants say so (#17).</summary>
        public static List<SpectralMatch> TruncationRowsToWrite(IEnumerable<SpectralMatch> psms, SearchParameters parameters) =>
            psms.Where(p => p != null
                    && (parameters.WriteDecoys || !p.IsDecoy)
                    && (parameters.WriteContaminants || !p.IsContaminant))
                .ToList();

        /// <summary>The Pass 3 chop acceptor for a precursor tolerance; its notch count drives the pooled FDR (#9, #15).</summary>
        private static MassDiffAcceptor GetChopAcceptor(MzLibUtil.Tolerance precursorMassTolerance) =>
            GetMassDiffAcceptor(precursorMassTolerance, TruncationChopAcceptorType, null);

        /// <summary>
        /// Dedups the search's PSMs to proteoform level (best PSM per FullSequence, mirroring AllProteoforms.psmtsv, #28),
        /// applies the permissive parent filter (#3), and expands each matched biopolymer into a parent — pipe-ambiguous
        /// proteoforms become separate parents (#2).
        /// </summary>
        private List<TruncationParent> BuildTruncationParents(IEnumerable<SpectralMatch> pass1Psms)
        {
            var parents = new List<TruncationParent>();

            // One PSM per proteoform, chosen the way FdrAnalysisEngine chooses its proteoform-level representative
            // (lowest PEP, then score), so the PSM kept is the one carrying the proteoform-level q-values. A
            // pipe-ambiguous match has a null FullSequence; keying it on its hypotheses keeps each ambiguous
            // proteoform separate rather than collapsing all of them into one null group (#2).
            IEnumerable<SpectralMatch> proteoforms = pass1Psms
                .Where(p => p != null)
                .OrderBy(p => p.FdrInfo?.PEP ?? double.MaxValue)
                .ThenByDescending(p => p)
                .GroupBy(ProteoformKey)
                .Select(g => g.First());

            foreach (SpectralMatch psm in proteoforms)
            {
                if (!PassesTruncationParentFilter(psm, TruncationParentQValueThreshold))
                {
                    continue;
                }

                // Record whether the PEP q-value (vs notch q-value) drove inclusion (perf logging).
                FdrInfo fdr = ParentFdrInfo(psm);
                if (fdr != null && fdr.PEP_QValue != 2)
                {
                    _pepQWasUsedForParents = true;
                }

                foreach (var hypothesis in psm.BestMatchingBioPolymersWithSetMods)
                {
                    if (hypothesis.SpecificBioPolymer is PeptideWithSetModifications proteoform)
                    {
                        parents.Add(new TruncationParent(proteoform, proteoform.Parent.Accession, psm, proteoform.Parent.IsDecoy));
                    }
                }
            }

            return parents;
        }

        /// <summary>Dedup key for a Pass 1 match: its FullSequence, or for a pipe-ambiguous match (null FullSequence) its sorted hypotheses.</summary>
        public static string ProteoformKey(SpectralMatch psm) =>
            psm.FullSequence ?? string.Join("|", psm.BestMatchingBioPolymersWithSetMods
                .Select(b => b.SpecificBioPolymer.FullSequence).Distinct().OrderBy(s => s, StringComparer.Ordinal));

        /// <summary>
        /// The FDR values the parent filter reads: proteoform-level (peptide) when FdrAnalysisEngine assigned them,
        /// else PSM-level. Proteoform-level values exist only on the engine's one representative per FullSequence;
        /// the other members, and every pipe-ambiguous match after the first (MetaMorpheus groups all null
        /// FullSequences together), keep the unassigned <see cref="FdrInfo()"/> sentinel of 2 in every q column.
        /// </summary>
        public static FdrInfo ParentFdrInfo(SpectralMatch psm)
        {
            FdrInfo peptide = psm.GetFdrInfo(peptideLevel: true);
            bool unassigned = peptide == null || (peptide.QValue == 2 && peptide.QValueNotch == 2 && peptide.PEP_QValue == 2);
            return unassigned ? psm.GetFdrInfo(peptideLevel: false) : peptide;
        }

        /// <summary>
        /// Permissive parent filter (#3): PEP q-value ≤ threshold when PEP was computed (sentinel 2 = not
        /// computed), otherwise notch q-value ≤ threshold. Uses proteoform-level (peptide) FDR, falling
        /// back to PSM-level where proteoform-level FDR was not assigned (see <see cref="ParentFdrInfo"/>).
        /// </summary>
        public static bool PassesTruncationParentFilter(SpectralMatch psm, double threshold)
        {
            FdrInfo fdr = ParentFdrInfo(psm);
            if (fdr == null)
            {
                return false;
            }

            return fdr.PEP_QValue != 2 ? fdr.PEP_QValue <= threshold : fdr.QValueNotch <= threshold;
        }
    }
}
