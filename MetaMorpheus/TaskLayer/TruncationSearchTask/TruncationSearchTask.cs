using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using EngineLayer.Truncation;
using MassSpectrometry;
using Omics;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    /// <summary>
    /// Identifies N- and C-terminally truncated proteoforms by re-searching MS2 scans against a
    /// fragment index built from the proteoforms found by an upstream top-down <see cref="SearchTask"/>.
    /// See docs/Truncation-Search.md for the full three-pass design (index + dual single-series scoring,
    /// terminus-directed chopping, pooled FDR/PEP).
    ///
    /// <see cref="RunSpecific"/> wires the Phase 1-3 engine pieces together: ingest the deduped Pass 1
    /// proteoforms (in-memory <see cref="TaskChainContext"/> hand-off, decision #1, with an
    /// AllProteoforms.psmtsv disk fallback), build target+decoy parents, run Pass 2 dual single-series
    /// scoring per file, chop the winners (Pass 3), pool with the inherited intact matches, run pooled
    /// FDR/PEP, and write AllTruncatedPSMs/AllTruncatedProteoforms.psmtsv.
    /// </summary>
    public class TruncationSearchTask : MetaMorpheusTask
    {
        public const string TruncatedPsmsFileName = "AllTruncatedPSMs.psmtsv";
        public const string TruncatedProteoformsFileName = "AllTruncatedProteoforms.psmtsv";

        public TruncationSearchTask() : base(MyTask.Truncation)
        {
            CommonParameters = new CommonParameters();
            TruncationSearchParameters = new TruncationSearchParameters();
        }

        public TruncationSearchParameters TruncationSearchParameters { get; set; }

        /// <summary>This task ingests the upstream search's PSM set from the task-chain context (#1).</summary>
        // [TomlIgnore] must be repeated on the override: Nett reads the attribute off the runtime type's own
        // property, so it does not inherit the one on MetaMorpheusTask and would otherwise serialize this
        // read-only property into the task TOML (and then fail to read it back).
        [Nett.TomlIgnore]
        public override bool ConsumesTaskChainContext => true;

        // Set during parent building; reported in the perf log (did PEP, not notch q, drive the #3 filter?).
        private bool _pepQWasUsedForParents;

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList,
            List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            MyTaskResults = new MyTaskResults(this);
            var wallStopwatch = System.Diagnostics.Stopwatch.StartNew();

            // The chop acceptor must map a precursor to a whole-residue chop at a notch (decisions #9, #15);
            // fail before any work if the configured type cannot (see ValidateChopAcceptor).
            MassDiffAcceptor chopAcceptor = GetChopAcceptor(CommonParameters.PrecursorMassTolerance);
            ValidateChopAcceptor(chopAcceptor, TruncationSearchParameters.MassDiffAcceptorType);

            // 1. Ingest the deduped Pass 1 proteoform-level matches (decision #1). The in-memory list
            //    keeps the matched proteoform's Protein + DigestionParams, which Pass 3 chopping and
            //    decoy generation need; the disk path reconstructs them.
            bool haveInMemoryProvenance = TryIngestFromContext(out List<SpectralMatch> pass1Psms);
            WarnOnUpstreamParameterMismatch();

            // Per-file Pass 1 PSM lookup, used only on the in-memory path to inherit intact matches as
            // full-length forms (#4a). Disk-ingested parents have no per-scan provenance.
            ILookup<string, SpectralMatch> pass1ByFile = haveInMemoryProvenance
                ? pass1Psms.Where(p => p != null).ToLookup(p => p.FullFilePath)
                : null;

            LoadModifications(taskId, out List<Modification> variableModifications,
                out List<Modification> fixedModifications, out List<string> localizableModificationTypes);

            var myFileManager = new MyFileManager(true);
            (CommonParameters FileParams, Ms2ScanWithSpecificMass[] Scans) LoadScans(int fileIndex)
            {
                string rawFilePath = currentRawFileList[fileIndex];
                CommonParameters fileParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[fileIndex]);
                MsDataFile dataFile = myFileManager.LoadFile(rawFilePath, fileParams);
                Ms2ScanWithSpecificMass[] scans = GetMs2Scans(dataFile, rawFilePath, fileParams)
                    .OrderBy(s => s.OneBasedScanNumber).ToArray();
                myFileManager.DoneWithFile(rawFilePath);
                return (fileParams, scans);
            }

            // The sequence-tag filter needs every file's scans before parents are built, and its per-scan
            // restriction is keyed on those scan objects, so on that path all files are loaded up front and
            // reused by the engine loop. Every other path loads one file at a time.
            bool tagFilter = TruncationSearchParameters.SeedParentsFromDatabase && TruncationSearchParameters.UseSequenceTagFilter;
            var preloaded = tagFilter ? Enumerable.Range(0, currentRawFileList.Count).Select(LoadScans).ToList() : null;

            // Build parents: database-seeded (optionally narrowed by the sequence-tag filter, decision: remove
            // the observed-parent requirement) takes precedence; otherwise seed from the upstream-identified
            // proteoforms (in-memory, else disk).
            List<TruncationParent> parents;
            Dictionary<Ms2ScanWithSpecificMass, HashSet<string>> perScanAllowed = null;
            if (TruncationSearchParameters.SeedParentsFromDatabase)
            {
                if (tagFilter)
                {
                    parents = BuildParentsFromDatabaseTagFiltered(taskId, dbFilenameList,
                        preloaded.SelectMany(f => f.Scans).ToList(), out perScanAllowed);
                }
                else
                {
                    parents = BuildParentsFromDatabase(dbFilenameList, variableModifications, fixedModifications, localizableModificationTypes);
                }
            }
            else
            {
                parents = haveInMemoryProvenance ? BuildParentsFromPsms(pass1Psms) : BuildParentsFromDisk();
            }

            // Per-scan tag restriction is applied in the engine only when both the tag filter and the per-scan
            // flag are on; otherwise null = no restriction (global-union or non-tag paths).
            if (TruncationSearchParameters.UsePerScanTagRestriction && perScanAllowed == null)
            {
                Warn("UsePerScanTagRestriction is set but has no effect: it requires both SeedParentsFromDatabase " +
                     "and UseSequenceTagFilter to be enabled. Proceeding without per-scan restriction.");
            }
            IReadOnlyDictionary<Ms2ScanWithSpecificMass, HashSet<string>> allowedByScan =
                TruncationSearchParameters.UsePerScanTagRestriction ? perScanAllowed : null;

            if (parents.Count == 0)
            {
                Warn("TruncationSearchTask found no proteoforms to seed the search; writing empty result files.");
                WriteOutputs(new List<SpectralMatch>(), OutputFolder, taskId);
                return MyTaskResults;
            }

            // 2. Augment with reverse decoys (decision #14).
            int targetAndPass1DecoyParents = parents.Count;
            parents = TruncationParentBuilder.AddReverseDecoys(parents);
            int decoysGenerated = parents.Count - targetAndPass1DecoyParents;

            var pooled = new List<SpectralMatch>();
            var internalPooled = new List<SpectralMatch>(); // 4b, filled per file
            // Diagnostic: winning parent's 0-based rank in each scan's index match-count ordering (validates
            // the engine's MaxCandidatesToScore cap). Written to CandidateRanks.tsv.
            var winnerRankRows = new List<string>();

            // Perf accumulators (docs/Truncation-Search.md, Benchmarking hook); parent/oversize counts are identical per file (shared list).
            var pass3Timings = new TruncationTimings();
            double pass2IndexSeconds = 0, pass2ScoringSeconds = 0;
            int totalMs2Scans = 0, indexedParents = 0, oversizeExcluded = 0;

            for (int i = 0; i < currentRawFileList.Count; i++)
            {
                // Load one file's MS2 scans at a time (unless the tag filter preloaded them), so only the file
                // being searched is held in memory.
                string rawFilePath = currentRawFileList[i];
                (CommonParameters fileParams, Ms2ScanWithSpecificMass[] scans) = preloaded != null ? preloaded[i] : LoadScans(i);
                if (preloaded != null)
                {
                    preloaded[i] = default; // release this file's scans once searched
                }

                // Every Pass 1 PSM on each scan, best first: a chimeric scan can carry several confident
                // intact matches, one per precursor, and the engine matches each precursor entry to its own (#4a).
                Dictionary<int, IReadOnlyList<SpectralMatch>> pass1PsmsByScan = pass1ByFile?[rawFilePath]
                    .GroupBy(p => p.ScanNumber)
                    .ToDictionary(g => g.Key, g => (IReadOnlyList<SpectralMatch>)g.OrderByDescending(p => p).ToList());

                // Pass 3 chop acceptor at this file's precursor tolerance, the same one Pass 2 uses (#9, #15).
                MassDiffAcceptor fileChopAcceptor = GetChopAcceptor(fileParams.PrecursorMassTolerance);

                var file = (Raw: rawFilePath, FileParams: fileParams, Scans: scans);
                var engine = new TruncationSearchEngine(parents, file.Scans, file.FileParams,
                    new TruncationAcceptor(file.FileParams.PrecursorMassTolerance), pass1PsmsByScan,
                    TruncationSearchParameters.MaxParentMass, allowedByScan);
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

                // Internal-fragment search on this file's scans (its FDR runs separately, in 4b). Collapsed per
                // file, since scan indices are only unique within a file.
                if (TruncationSearchParameters.SearchInternalTruncations)
                {
                    internalPooled.AddRange(TruncationPass3.CollapseDuplicateTruncations(InternalTruncationSearch.Run(
                        parents, file.Scans, file.FileParams, fileChopAcceptor,
                        TruncationSearchParameters.InternalMinIonsPerTerminus, TruncationSearchParameters.MaxParentMass)));
                }
            }

            // 4. Pooled FDR + PEP over (full-length + truncation) PSMs (decisions #15, #18). An empty pool
            //    (no winners and no inherited intact matches) skips FDR and writes header-only files.
            var fdrStopwatch = System.Diagnostics.Stopwatch.StartNew();
            List<SpectralMatch> withFdr = pooled.Count == 0
                ? pooled
                : TruncationFdr.RunPooledFdr(pooled, CommonParameters,
                    chopAcceptor.NumNotches, FileSpecificParameters, taskId, OutputFolder);
            double fdrPepSeconds = fdrStopwatch.Elapsed.TotalSeconds;

            // 4b. Internal-fragment search (separate per-class FDR): a real internal fragment shares no ladder
            //     with its parent, so it needs a direct search scored on its own b/y ladder (propensity-weighted)
            //     and judged among internal targets/decoys only — not pooled against the easy intact threshold.
            if (TruncationSearchParameters.SearchInternalTruncations)
            {
                List<SpectralMatch> internalWithFdr = internalPooled.Count == 0
                    ? internalPooled
                    : TruncationFdr.RunPooledFdr(internalPooled, CommonParameters,
                        chopAcceptor.NumNotches, FileSpecificParameters, taskId, OutputFolder);

                string internalPath = Path.Combine(OutputFolder, "AllInternalTruncations.psmtsv");
                TruncationOutput.WritePsms(internalWithFdr, internalPath);
                FinishedWritingFile(internalPath, new List<string> { taskId });
            }

            // 5. Write the two result files (decisions #16, #17). No q/PEP cutoff; decoys/contaminants included.
            WriteOutputs(withFdr, OutputFolder, taskId);

            // Diagnostic side file: per-scan winning-parent rank in the index match-count ordering.
            // Opt-in (gated on the same PerfLogPath that enables benchmarking), registered via
            // FinishedWritingFile, and best-effort so a diagnostic-write failure can't abort an otherwise
            // successful task (#14).
            if (!string.IsNullOrWhiteSpace(TruncationSearchParameters.PerfLogPath))
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
            }

            // 6. Optional perf-log row (docs/Truncation-Search.md, Benchmarking hook). No-op unless a path is configured.
            if (!string.IsNullOrWhiteSpace(TruncationSearchParameters.PerfLogPath))
            {
                AppendPerfLog(withFdr, OutputFolder, currentRawFileList.Count, totalMs2Scans,
                    wallStopwatch.Elapsed.TotalSeconds, pass2IndexSeconds, pass2ScoringSeconds, pass3Timings,
                    fdrPepSeconds, indexedParents, oversizeExcluded, decoysGenerated);
            }

            return MyTaskResults;
        }

        /// <summary>Computes the per-run metrics from the written set and appends one perf_log.tsv row.</summary>
        private void AppendPerfLog(List<SpectralMatch> withFdr, string outputFolder, int nRawFiles, int totalMs2Scans,
            double wallSeconds, double pass2IndexSeconds, double pass2ScoringSeconds, TruncationTimings pass3Timings,
            double fdrPepSeconds, int indexedParents, int oversizeExcluded, int decoysGenerated)
        {
            List<SpectralMatch> rows = withFdr.Where(p => p != null).ToList();
            int Truncations(string kind) => rows.Count(p => DescriptionOf(p).Contains(kind));
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
                PerfLogger.Append(TruncationSearchParameters.PerfLogPath, metrics);
            }
            catch (System.Exception ex)
            {
                Warn("Failed to append perf-log row: " + ex.Message);
            }
        }

        private static string DescriptionOf(SpectralMatch psm)
        {
            var pwsm = psm.BestMatchingBioPolymersWithSetMods.FirstOrDefault()?.SpecificBioPolymer as PeptideWithSetModifications;
            return pwsm?.Description ?? string.Empty;
        }

        /// <summary>
        /// Writes AllTruncatedPSMs.psmtsv and AllTruncatedProteoforms.psmtsv (always, even when empty). No q/PEP
        /// cutoff; decoy and contaminant rows follow WriteDecoys/WriteContaminants (#17), as in PostSearchAnalysisTask.
        /// </summary>
        private void WriteOutputs(List<SpectralMatch> psms, string outputFolder, string taskId)
        {
            List<SpectralMatch> toWrite = RowsToWrite(psms, TruncationSearchParameters);

            string psmsPath = Path.Combine(outputFolder, TruncatedPsmsFileName);
            string proteoformsPath = Path.Combine(outputFolder, TruncatedProteoformsFileName);
            TruncationOutput.WritePsms(toWrite, psmsPath, TruncationSearchParameters.ModsToWriteSelection);
            TruncationOutput.WriteProteoforms(toWrite, proteoformsPath, TruncationSearchParameters.ModsToWriteSelection);
            FinishedWritingFile(psmsPath, new List<string> { taskId });
            FinishedWritingFile(proteoformsPath, new List<string> { taskId });
        }

        /// <summary>The PSMs that reach the output files: decoys and contaminants only when WriteDecoys / WriteContaminants say so (#17).</summary>
        public static List<SpectralMatch> RowsToWrite(IEnumerable<SpectralMatch> psms, TruncationSearchParameters parameters) =>
            psms.Where(p => p != null
                    && (parameters.WriteDecoys || !p.IsDecoy)
                    && (parameters.WriteContaminants || !p.IsContaminant))
                .ToList();

        /// <summary>The Pass 3 chop acceptor for a precursor tolerance; its notch count drives the pooled FDR (#9, #15).</summary>
        private MassDiffAcceptor GetChopAcceptor(MzLibUtil.Tolerance precursorMassTolerance) =>
            SearchTask.GetMassDiffAcceptor(precursorMassTolerance,
                TruncationSearchParameters.MassDiffAcceptorType, TruncationSearchParameters.CustomMdac);

        /// <summary>
        /// Chopping only works with an acceptor that accepts the observed monoisotopic precursor at a fixed set of
        /// notches: an exact window or the 1/2/3 mm and ±3 mm notch sets (a Custom string must parse to one of
        /// these). Open and interval acceptors (Open, ModOpen) accept a mass shift, so the chopper would stop early
        /// and absorb it (#9 rules that out); the MostAbundant types expect the most-abundant mass, not the
        /// monoisotopic precursor the chopper passes, and never match.
        /// </summary>
        public static void ValidateChopAcceptor(MassDiffAcceptor chopAcceptor, MassDiffAcceptorType type)
        {
            if (chopAcceptor is DotMassDiffAcceptor or SinglePpmAroundZeroSearchMode or SingleAbsoluteAroundZeroSearchMode)
            {
                return;
            }

            throw new MetaMorpheusException($"TruncationSearchTask cannot chop with MassDiffAcceptorType {type} " +
                $"({chopAcceptor.GetType().Name}). Use Exact, OneMM, TwoMM, ThreeMM, PlusOrMinusThreeMM, or a Custom dot or ppm/da-around-zero acceptor.");
        }

        /// <summary>
        /// The task searches with its own CommonParameters, which default to bottom-up settings. Warns when the
        /// settings that decide MS2 deconvolution, precursor matching and fragment matching differ from the upstream
        /// search whose matches it consumes (a top-down run usually wants them identical).
        /// </summary>
        private void WarnOnUpstreamParameterMismatch()
        {
            if (TaskChainContext == null)
            {
                return;
            }

            string upstreamId = TruncationSearchParameters.UpstreamSearchTaskId;
            bool found = upstreamId != null
                ? TaskChainContext.TryGet(TaskChainContext.CommonParametersKey(upstreamId), out CommonParameters upstream)
                : TaskChainContext.TryGetMostRecent(out upstream);
            if (!found || upstream == null)
            {
                return;
            }

            List<string> differences = DescribeParameterDifferences(upstream, CommonParameters);
            if (differences.Count > 0)
            {
                Warn("TruncationSearchTask settings differ from the upstream search it consumes: " + string.Join("; ", differences)
                    + ". Copy the search's CommonParameters into the truncation task unless this is intended.");
            }
        }

        /// <summary>Human-readable "name: upstream vs this" differences in the settings that decide deconvolution and matching.</summary>
        public static List<string> DescribeParameterDifferences(CommonParameters upstream, CommonParameters mine)
        {
            var settings = new (string Name, Func<CommonParameters, object> Get)[]
            {
                ("PrecursorMassTolerance", c => c.PrecursorMassTolerance),
                ("ProductMassTolerance", c => c.ProductMassTolerance),
                ("DoPrecursorDeconvolution", c => c.DoPrecursorDeconvolution),
                ("UseProvidedPrecursorInfo", c => c.UseProvidedPrecursorInfo),
                ("PrecursorDeconvolution max charge", c => c.PrecursorDeconvolutionParameters?.MaxAssumedChargeState),
                ("ProductDeconvolution max charge", c => c.ProductDeconvolutionParameters?.MaxAssumedChargeState),
                ("DissociationType", c => c.DissociationType),
            };

            var differences = new List<string>();
            foreach ((string name, Func<CommonParameters, object> get) in settings)
            {
                string theirs = get(upstream)?.ToString() ?? "null";
                string ours = get(mine)?.ToString() ?? "null";
                if (theirs != ours)
                {
                    differences.Add($"{name} {ours} (search: {theirs})");
                }
            }

            return differences;
        }

        /// <summary>
        /// Primary ingest (decision #1): pull the upstream SearchTask's resolved, FDR'd PSM set from the
        /// shared in-memory <see cref="TaskChainContext"/>. These carry the matched Protein/DigestionParams
        /// that Pass 3 chopping and decoy generation require.
        /// </summary>
        private bool TryIngestFromContext(out List<SpectralMatch> pass1Psms)
        {
            pass1Psms = null;
            if (TaskChainContext == null)
            {
                return false;
            }

            string upstreamId = TruncationSearchParameters.UpstreamSearchTaskId;

            // Explicit upstream id when configured; otherwise fall back to the most recent upstream
            // result in the run list (so a [Search, Truncation] CMD run needs no exact-id wiring).
            bool found = upstreamId != null
                ? TaskChainContext.TryGet(upstreamId, out pass1Psms)
                : TaskChainContext.TryGetMostRecent(out pass1Psms);

            return found && pass1Psms != null;
        }

        /// <summary>
        /// Dedups the Pass 1 PSMs to proteoform level (best PSM per FullSequence, mirroring
        /// AllProteoforms.psmtsv, #28), applies the permissive parent filter (#3), and expands each
        /// matched biopolymer into a parent — pipe-ambiguous proteoforms become separate parents (#2).
        /// </summary>
        private List<TruncationParent> BuildParentsFromPsms(IEnumerable<SpectralMatch> pass1Psms)
        {
            double threshold = TruncationSearchParameters.ParentQValueThreshold;
            var parents = new List<TruncationParent>();

            IEnumerable<SpectralMatch> proteoforms = pass1Psms
                .Where(p => p != null)
                .OrderByDescending(p => p)
                .GroupBy(p => p.FullSequence)
                .Select(g => g.First());

            foreach (SpectralMatch psm in proteoforms)
            {
                if (!PassesParentFilter(psm, threshold))
                {
                    continue;
                }

                // Record whether the PEP q-value (vs notch q-value) drove inclusion (perf logging).
                // By design this reflects the in-memory parent path only; the disk/db-seeded paths leave the
                // metric false because PsmFromTsv's PEP column is ambiguous there (#15, intentional).
                FdrInfo fdr = psm.GetFdrInfo(peptideLevel: true) ?? psm.GetFdrInfo(peptideLevel: false);
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

        /// <summary>
        /// Permissive parent filter (#3): PEP q-value ≤ threshold when PEP was computed (sentinel 2 = not
        /// computed), otherwise notch q-value ≤ threshold. Uses proteoform-level (peptide) FDR, falling
        /// back to PSM-level if peptide-level FDR was not assigned.
        /// </summary>
        private static bool PassesParentFilter(SpectralMatch psm, double threshold)
        {
            FdrInfo fdr = psm.GetFdrInfo(peptideLevel: true) ?? psm.GetFdrInfo(peptideLevel: false);
            if (fdr == null)
            {
                return false;
            }

            return fdr.PEP_QValue != 2 ? fdr.PEP_QValue <= threshold : fdr.QValueNotch <= threshold;
        }

        /// <summary>
        /// Database-seeded parents: generate theoretical full-length proteoforms from the protein database the
        /// same way a top-down search does — load each protein (targets only; the truncation search makes its
        /// own reverse decoys) and digest it with the configured fixed/variable mods, so XML-annotated PTMs and
        /// variable mods ride along. This removes the requirement that a truncation's parent be identified
        /// intact (the dominant recall ceiling). Oversized parents are filtered later by the engine's
        /// <see cref="TruncationSearchParameters.MaxParentMass"/>.
        /// </summary>
        private List<TruncationParent> BuildParentsFromDatabase(List<DbForTask> dbFilenameList, List<Modification> variableModifications,
            List<Modification> fixedModifications, List<string> localizableModificationTypes)
        {
            var parents = new List<TruncationParent>();
            foreach (DbForTask db in dbFilenameList)
            {
                IEnumerable<Protein> proteins = DatabaseLoadingEngine.LoadProteinDb(
                    db.FilePath, generateTargets: true, decoyType: DecoyType.None,
                    localizableModificationTypes, db.IsContaminant, out _, out _, CommonParameters);

                foreach (Protein protein in proteins)
                {
                    if (protein.IsDecoy)
                    {
                        continue;
                    }

                    foreach (PeptideWithSetModifications proteoform in protein
                        .Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications)
                        .OfType<PeptideWithSetModifications>())
                    {
                        parents.Add(new TruncationParent(proteoform, protein.Accession, protein.Accession, isDecoy: false));
                    }
                }
            }

            return parents;
        }

        /// <summary>
        /// Database-seeded parents, narrowed by the sequence-tag filter (docs/Truncation-Search.md, Sequence-tag filtering): build a k-mer index
        /// over the database, extract de-novo tags from every scan, take the UNION of tag-supported proteins
        /// (global-union variant), and digest only those into theoretical parents — keeping the parent set
        /// far smaller than the whole database without restricting candidates per scan. Tag extraction (over
        /// scans) and digestion (over candidate proteins) both run in parallel.
        /// </summary>
        private List<TruncationParent> BuildParentsFromDatabaseTagFiltered(string taskId, List<DbForTask> dbFilenameList,
            IReadOnlyList<Ms2ScanWithSpecificMass> allScans,
            out Dictionary<Ms2ScanWithSpecificMass, HashSet<string>> perScanAllowedAccessions)
        {
            LoadModifications(taskId, out List<Modification> variableModifications,
                out List<Modification> fixedModifications, out List<string> localizableModificationTypes);

            var proteins = new List<Protein>();
            foreach (DbForTask db in dbFilenameList)
            {
                proteins.AddRange(DatabaseLoadingEngine.LoadProteinDb(
                        db.FilePath, generateTargets: true, decoyType: DecoyType.None,
                        localizableModificationTypes, db.IsContaminant, out _, out _, CommonParameters)
                    .Where(p => !p.IsDecoy));
            }

            int tagLength = TruncationSearchParameters.TagLength;
            int minTagHits = TruncationSearchParameters.MinTagHits;

            // TagLength and MinTagHits are user-settable in the TOML and reach ProteinTagIndex unchecked.
            // A tag length below 1 throws ArgumentOutOfRangeException from its constructor -- a CMD stack
            // trace today, and the crash-report dialog the moment this task is wired into the GUI. An
            // impossible configuration gets a Warn and the default here, which is the house channel for it.
            if (tagLength < 1)
            {
                Warn($"TagLength must be at least 1; {tagLength} was given. Using the default of " +
                     $"{TaskLayer.TruncationSearchParameters.DefaultTagLength} instead.");
                tagLength = TaskLayer.TruncationSearchParameters.DefaultTagLength;
            }
            if (minTagHits < 1)
            {
                Warn($"MinTagHits must be at least 1; {minTagHits} was given. Using the default of " +
                     $"{TaskLayer.TruncationSearchParameters.DefaultMinTagHits} instead -- a value below 1 asks no protein to carry any tag, " +
                     $"which selects the whole database and defeats the filter.");
                minTagHits = TaskLayer.TruncationSearchParameters.DefaultMinTagHits;
            }

            int threads = Math.Max(1, CommonParameters.MaxThreadsToUsePerFile);
            var index = new ProteinTagIndex(proteins, tagLength, threads);
            MzLibUtil.Tolerance productTolerance = CommonParameters.ProductMassTolerance;
            var parallelOptions = new ParallelOptions { MaxDegreeOfParallelism = threads };
            const int maxFragmentPeaksForTags = 120; // tags come from the most intense peaks; noise peaks excluded

            // Per scan: the proteins its de-novo tags support (recorded for the per-scan restriction); the union
            // is what we digest into parents so the shared index covers every scan's candidates.
            var perThreadCandidates = new ConcurrentBag<HashSet<int>>();
            var perScanAccessions = new ConcurrentDictionary<Ms2ScanWithSpecificMass, HashSet<string>>();
            Parallel.ForEach(allScans, parallelOptions, () => new HashSet<int>(), (scan, _, local) =>
            {
                if (scan.ExperimentalFragments != null && scan.ExperimentalFragments.Length > 0)
                {
                    // Only the most intense fragment peaks feed the tag extractor — low-intensity noise peaks
                    // create spurious residue-sized gaps that otherwise let tags match nearly every protein.
                    List<double> masses = scan.ExperimentalFragments
                        .OrderByDescending(e => e.TotalIntensity)
                        .Take(maxFragmentPeaksForTags)
                        .Select(e => e.MonoisotopicMass)
                        .ToList();
                    HashSet<string> tags = SequenceTagExtractor.ExtractTags(masses, productTolerance, tagLength);
                    var scanAccessions = new HashSet<string>();
                    foreach (int id in index.GetCandidateProteinIds(tags, minTagHits))
                    {
                        local.Add(id);
                        scanAccessions.Add(proteins[id].Accession);
                    }
                    perScanAccessions[scan] = scanAccessions;
                }
                return local;
            }, local => perThreadCandidates.Add(local));

            perScanAllowedAccessions = new Dictionary<Ms2ScanWithSpecificMass, HashSet<string>>(perScanAccessions);

            var candidateProteinIds = new HashSet<int>();
            foreach (HashSet<int> set in perThreadCandidates) candidateProteinIds.UnionWith(set);

            Warn($"Sequence-tag filter selected {candidateProteinIds.Count} of {proteins.Count} database proteins " +
                 $"(tag length {tagLength}, min {minTagHits} tag hits).");

            // Digest only the candidate proteins (in parallel) into theoretical full-length parents.
            //
            // One slot per candidate protein, concatenated in protein-id order, rather than a
            // ConcurrentBag of per-thread lists. This is the first parent-building path that is
            // parallel -- BuildParentsFromDatabase is a serial nested foreach, so its parent order is a
            // function of the database -- and parent order is observable: TruncationSearchEngine sorts
            // parents by mass with OrderBy, which is stable, so equal-mass parents keep their arrival
            // order, and Pass 2 breaks a score tie with `nScore > bestScore`, giving the win to whichever
            // arrived first. Left as a bag, two runs over the same database could name different
            // proteins for the same scan.
            int[] orderedCandidateIds = candidateProteinIds.OrderBy(id => id).ToArray();
            var parentsByCandidate = new List<TruncationParent>[orderedCandidateIds.Length];
            Parallel.For(0, orderedCandidateIds.Length, parallelOptions, i =>
            {
                Protein protein = proteins[orderedCandidateIds[i]];
                var local = new List<TruncationParent>();
                foreach (PeptideWithSetModifications proteoform in protein
                    .Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications)
                    .OfType<PeptideWithSetModifications>())
                {
                    local.Add(new TruncationParent(proteoform, protein.Accession, protein.Accession, isDecoy: false));
                }
                parentsByCandidate[i] = local;
            });

            var parents = new List<TruncationParent>();
            foreach (List<TruncationParent> list in parentsByCandidate) parents.AddRange(list);
            return parents;
        }

        /// <summary>
        /// Disk fallback (decision #1): build parents from an AllProteoforms.psmtsv. Disk rows carry no
        /// Protein, so a synthetic single-chain Protein is reconstructed per parent so Pass 3 chopping and
        /// reverse-decoy generation still work. The proteoform is placed at its real start in the protein
        /// (from the row's "Start and End Residues In Protein", with unknown residues before it), so truncation
        /// coordinates stay protein-relative (#13) as on the in-memory path; mods are parsed from the
        /// FullSequence. Pipe-ambiguous rows expand into separate parents (#2). Intact-match inheritance
        /// (#4a) is unavailable on this path (no per-scan provenance).
        /// </summary>
        private List<TruncationParent> BuildParentsFromDisk()
        {
            string path = TruncationSearchParameters.Pass1ProteoformsFilePath;
            if (path == null || !File.Exists(path))
            {
                return new List<TruncationParent>();
            }

            List<PsmFromTsv> rows = SpectrumMatchTsvReader.ReadPsmTsv(path, out _);
            double threshold = TruncationSearchParameters.ParentQValueThreshold;
            var digestionParams = CommonParameters.DigestionParams as DigestionParams;
            var parents = new List<TruncationParent>();

            foreach (PsmFromTsv row in rows)
            {
                if (!PassesDiskParentFilter(row, threshold))
                {
                    continue;
                }

                string[] sequences = row.FullSequence.Split('|');
                string[] accessions = (row.Accession ?? string.Empty).Split('|');
                List<int> starts = ParseStartResidues(row.StartAndEndResiduesInProtein);

                for (int i = 0; i < sequences.Length; i++)
                {
                    string fullSequence = sequences[i];
                    string accession = i < accessions.Length && accessions[i].Length > 0 ? accessions[i]
                        : accessions.Length > 0 && accessions[0].Length > 0 ? accessions[0] : "UNKNOWN";
                    int start = i < starts.Count ? starts[i] : starts.Count > 0 ? starts[0] : 1;

                    parents.Add(new TruncationParent(BuildDiskProteoform(fullSequence, accession, start, row.IsDecoy, digestionParams),
                        accession, row, row.IsDecoy));
                }
            }

            return parents;
        }

        /// <summary>
        /// A disk parent on a synthetic protein: unknown residues ('X') up to <paramref name="oneBasedStart"/>, then the
        /// proteoform, so its coordinates (and every truncation chopped from it) are the protein's own (#13).
        /// </summary>
        public static PeptideWithSetModifications BuildDiskProteoform(string fullSequence, string accession, int oneBasedStart,
            bool isDecoy, DigestionParams digestionParams)
        {
            string baseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSequence);
            var protein = new Protein(new string('X', oneBasedStart - 1) + baseSequence, accession, isDecoy: isDecoy);
            return new PeptideWithSetModifications(fullSequence, GlobalVariables.AllModsKnownDictionary,
                digestionParams: digestionParams, p: protein,
                oneBasedStartResidueInProtein: oneBasedStart, oneBasedEndResidueInProtein: oneBasedStart + baseSequence.Length - 1,
                cleavageSpecificity: Omics.Digestion.CleavageSpecificity.Full);
        }

        /// <summary>One-based start residues from a "[2 to 120]|[5 to 123]" cell, one per pipe-separated alternative.</summary>
        public static List<int> ParseStartResidues(string startAndEndResidues) =>
            System.Text.RegularExpressions.Regex.Matches(startAndEndResidues ?? string.Empty, @"\[(\d+) to \d+\]")
                .Select(m => int.Parse(m.Groups[1].Value))
                .Where(start => start >= 1)
                .ToList();

        /// <summary>Permissive parent filter (#3) for disk rows (PsmFromTsv stores PEP_QValue/QValueNotch).</summary>
        private static bool PassesDiskParentFilter(PsmFromTsv row, double threshold)
        {
            // 0 is treated as "PEP not computed" alongside the 2 sentinel ON PURPOSE for disk rows:
            // PsmFromTsv defaults an absent/empty PEP_QValue column to 0, so a disk 0 is ambiguous (absent
            // column vs a genuine best-possible 0.0) and the two cannot be told apart here. Falling back to
            // the notch q-value when PEP==0 is the safe choice. Do NOT "simplify" this to `!= 2` only -- that
            // would pass absent-PEP rows as if q=0 (reviewed: findings #11/#15, intentional wontfix).
            if (row.PEP_QValue != 0 && row.PEP_QValue != 2)
            {
                return row.PEP_QValue <= threshold;
            }

            return (row.QValueNotch ?? row.QValue) <= threshold;
        }
    }
}
