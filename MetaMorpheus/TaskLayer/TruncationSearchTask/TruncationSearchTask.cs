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
    /// See 01_Architecture.md for the full three-pass design (index + dual single-series scoring,
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

        // Set during parent building; reported in the perf log (did PEP, not notch q, drive the #3 filter?).
        private bool _pepQWasUsedForParents;

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList,
            List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            MyTaskResults = new MyTaskResults(this);
            var wallStopwatch = System.Diagnostics.Stopwatch.StartNew();

            // 1. Ingest the deduped Pass 1 proteoform-level matches (decision #1). The in-memory list
            //    keeps the matched proteoform's Protein + DigestionParams, which Pass 3 chopping and
            //    decoy generation need; the disk path reconstructs them.
            bool haveInMemoryProvenance = TryIngestFromContext(out List<SpectralMatch> pass1Psms);

            // Per-file Pass 1 PSM lookup, used only on the in-memory path to inherit intact matches as
            // full-length forms (#4a). Disk-ingested parents have no per-scan provenance.
            ILookup<string, SpectralMatch> pass1ByFile = haveInMemoryProvenance
                ? pass1Psms.Where(p => p != null).ToLookup(p => p.FullFilePath)
                : null;

            // Load every file's MS2 scans up front (kept in memory) plus its intact-skip dicts. Hoisted so the
            // sequence-tag filter can see all scans before parents are built; the engine loop reuses them.
            var myFileManager = new MyFileManager(true);
            var loaded = new List<(string Raw, CommonParameters FileParams, Ms2ScanWithSpecificMass[] Scans,
                Dictionary<int, SpectralMatch> Pass1PsmByScan, Dictionary<int, double> Pass1MassByScan)>();
            for (int i = 0; i < currentRawFileList.Count; i++)
            {
                string rawFilePath = currentRawFileList[i];
                CommonParameters fileParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[i]);

                MsDataFile dataFile = myFileManager.LoadFile(rawFilePath, fileParams);
                Ms2ScanWithSpecificMass[] scans = GetMs2Scans(dataFile, rawFilePath, fileParams)
                    .OrderBy(s => s.OneBasedScanNumber).ToArray();
                myFileManager.DoneWithFile(rawFilePath);

                Dictionary<int, SpectralMatch> pass1PsmByScan = pass1ByFile?[rawFilePath]
                    .GroupBy(p => p.ScanNumber)
                    .ToDictionary(g => g.Key, g => g.OrderByDescending(p => p).First());
                Dictionary<int, double> pass1MassByScan = pass1PsmByScan?
                    .Where(kv => kv.Value.BioPolymerWithSetModsMonoisotopicMass.HasValue)
                    .ToDictionary(kv => kv.Key, kv => kv.Value.BioPolymerWithSetModsMonoisotopicMass.Value);

                loaded.Add((rawFilePath, fileParams, scans, pass1PsmByScan, pass1MassByScan));
            }

            // Build parents: database-seeded (optionally narrowed by the sequence-tag filter, decision: remove
            // the observed-parent requirement) takes precedence; otherwise seed from the upstream-identified
            // proteoforms (in-memory, else disk).
            List<TruncationParent> parents;
            Dictionary<Ms2ScanWithSpecificMass, HashSet<string>> perScanAllowed = null;
            if (TruncationSearchParameters.SeedParentsFromDatabase)
            {
                if (TruncationSearchParameters.UseSequenceTagFilter)
                {
                    parents = BuildParentsFromDatabaseTagFiltered(taskId, dbFilenameList,
                        loaded.SelectMany(f => f.Scans).ToList(), out perScanAllowed);
                }
                else
                {
                    parents = BuildParentsFromDatabase(taskId, dbFilenameList);
                }
            }
            else
            {
                parents = haveInMemoryProvenance ? BuildParentsFromPsms(pass1Psms) : BuildParentsFromDisk();
            }

            // Per-scan tag restriction is applied in the engine only when both the tag filter and the per-scan
            // flag are on; otherwise null = no restriction (global-union or non-tag paths).
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

            // 3. Notch acceptor for Pass 3 chopping (decision #80); its notch count drives the pooled FDR.
            MassDiffAcceptor chopAcceptor = SearchTask.GetMassDiffAcceptor(
                CommonParameters.PrecursorMassTolerance,
                TruncationSearchParameters.MassDiffAcceptorType,
                TruncationSearchParameters.CustomMdac);

            var pooled = new List<SpectralMatch>();
            // Diagnostic: winning parent's 0-based rank in each scan's index match-count ordering (validates
            // the engine's MaxCandidatesToScore cap). Written to CandidateRanks.tsv.
            var winnerRankRows = new List<string>();

            // Perf accumulators (03_Benchmarks); parent/oversize counts are identical per file (shared list).
            var pass3Timings = new TruncationTimings();
            double pass2IndexSeconds = 0, pass2ScoringSeconds = 0;
            int totalMs2Scans = 0, indexedParents = 0, oversizeExcluded = 0;

            foreach (var file in loaded)
            {
                var engine = new TruncationSearchEngine(parents, file.Scans, file.FileParams,
                    new TruncationAcceptor(file.FileParams.PrecursorMassTolerance), file.Pass1MassByScan,
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
                        TruncationPsm psm = TruncationPass3.ScoreTruncation(selection, file.FileParams, chopAcceptor, pass3Timings);
                        if (psm != null)
                        {
                            truncationPsms.Add(psm);
                        }
                        winnerRankRows.Add(string.Join("\t", Path.GetFileNameWithoutExtension(file.Raw),
                            selection.Scan.OneBasedScanNumber, selection.CandidateRank, selection.CandidatePoolSize,
                            selection.Score.ToString("F3", System.Globalization.CultureInfo.InvariantCulture),
                            psm != null ? 1 : 0));
                    }
                    else if (selection.Outcome == TruncationScanOutcome.IntactInherited
                             && file.Pass1PsmByScan != null
                             && file.Pass1PsmByScan.TryGetValue(selection.Scan.OneBasedScanNumber, out SpectralMatch intactHit))
                    {
                        pooled.Add(TruncationPass3.InheritAsFullLength(intactHit, selection.Scan, selection.ScanIndex, file.FileParams));
                    }
                }

                // Collapse duplicate truncations (same truncated FullSequence + scan) before pooling (#10).
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

            // 5. Write the two result files (decisions #16, #17). No q/PEP cutoff; decoys/contaminants included.
            WriteOutputs(withFdr, OutputFolder, taskId);

            // Diagnostic side file: per-scan winning-parent rank in the index match-count ordering.
            string ranksPath = Path.Combine(OutputFolder, "CandidateRanks.tsv");
            File.WriteAllLines(ranksPath,
                new[] { "File\tScanNumber\tCandidateRank\tCandidatePoolSize\tScore\tProducedTruncation" }.Concat(winnerRankRows));

            // 6. Optional perf-log row (03_Benchmarks). No-op unless a path is configured.
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

        /// <summary>Writes AllTruncatedPSMs.psmtsv and AllTruncatedProteoforms.psmtsv (always, even when empty).</summary>
        private void WriteOutputs(List<SpectralMatch> psms, string outputFolder, string taskId)
        {
            string psmsPath = Path.Combine(outputFolder, TruncatedPsmsFileName);
            string proteoformsPath = Path.Combine(outputFolder, TruncatedProteoformsFileName);
            TruncationOutput.WritePsms(psms, psmsPath);
            TruncationOutput.WriteProteoforms(psms, proteoformsPath);
            FinishedWritingFile(psmsPath, new List<string> { taskId });
            FinishedWritingFile(proteoformsPath, new List<string> { taskId });
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
        private List<TruncationParent> BuildParentsFromDatabase(string taskId, List<DbForTask> dbFilenameList)
        {
            LoadModifications(taskId, out List<Modification> variableModifications,
                out List<Modification> fixedModifications, out List<string> localizableModificationTypes);

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
        /// Database-seeded parents, narrowed by the sequence-tag filter (doc §11.2.4): build a k-mer index
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
            var perThreadParents = new ConcurrentBag<List<TruncationParent>>();
            Parallel.ForEach(candidateProteinIds, parallelOptions, () => new List<TruncationParent>(), (id, _, local) =>
            {
                Protein protein = proteins[id];
                foreach (PeptideWithSetModifications proteoform in protein
                    .Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications)
                    .OfType<PeptideWithSetModifications>())
                {
                    local.Add(new TruncationParent(proteoform, protein.Accession, protein.Accession, isDecoy: false));
                }
                return local;
            }, local => perThreadParents.Add(local));

            var parents = new List<TruncationParent>();
            foreach (List<TruncationParent> list in perThreadParents) parents.AddRange(list);
            return parents;
        }

        /// <summary>
        /// Disk fallback (decision #1): build parents from an AllProteoforms.psmtsv. Disk rows carry no
        /// Protein, so a synthetic single-chain Protein is reconstructed per parent (whole sequence,
        /// start=1) so Pass 3 chopping and reverse-decoy generation still work; mods are parsed from the
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

                for (int i = 0; i < sequences.Length; i++)
                {
                    string fullSequence = sequences[i];
                    string accession = i < accessions.Length && accessions[i].Length > 0 ? accessions[i]
                        : accessions.Length > 0 && accessions[0].Length > 0 ? accessions[0] : "UNKNOWN";
                    string baseSequence = IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSequence);

                    var protein = new Protein(baseSequence, accession, isDecoy: row.IsDecoy);
                    var proteoform = new PeptideWithSetModifications(fullSequence, GlobalVariables.AllModsKnownDictionary,
                        digestionParams: digestionParams, p: protein,
                        oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: baseSequence.Length,
                        cleavageSpecificity: Omics.Digestion.CleavageSpecificity.Full);

                    parents.Add(new TruncationParent(proteoform, accession, row, row.IsDecoy));
                }
            }

            return parents;
        }

        /// <summary>Permissive parent filter (#3) for disk rows (PsmFromTsv stores PEP_QValue/QValueNotch).</summary>
        private static bool PassesDiskParentFilter(PsmFromTsv row, double threshold)
        {
            if (row.PEP_QValue != 0 && row.PEP_QValue != 2)
            {
                return row.PEP_QValue <= threshold;
            }

            return (row.QValueNotch ?? row.QValue) <= threshold;
        }
    }
}
