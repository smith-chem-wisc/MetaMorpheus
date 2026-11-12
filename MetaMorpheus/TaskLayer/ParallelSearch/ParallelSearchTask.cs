#nullable enable
using Chromatography.RetentionTimeCalibration;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using EngineLayer.Indexing;
using EngineLayer.ParallelSearch;
using EngineLayer.ParallelSearch.FdrAlignment;
using EngineLayer.SpectrumMatch;
using EngineLayer.Util;
using MassSpectrometry;
using Nett;
using Omics;
using Omics.Modifications;
using Omics.SpectrumMatch;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Channels;
using System.Threading.Tasks;
using TaskLayer.ParallelSearch.Analysis;
using TaskLayer.ParallelSearch.Analysis.Collectors;
using TaskLayer.ParallelSearch.Util;

namespace TaskLayer.ParallelSearch;
public class ParallelSearchTask : SearchTask
{
    private readonly object _progressLock = new object();
    private readonly object _countLock = new object();
    private TransientDatabaseResultsManager? _resultsManager;
    private Channel<(TransientDatabaseContext Context, TransientDatabaseMetrics Metrics)>? _completedDatabaseWriteChannel;
    private Task? _completedDatabaseWriterTask;
    private int _cachedAtStart;
    private int _initialFinishedCount;
    private int _completedThisRun;

    public ParallelSearchTask() : base(MyTask.ParallelSearch)
    {
        // Initialize with appropriate defaults
        SearchParameters = new ParallelSearchParameters();
        CommonParameters = new(taskDescriptor: "ParallelSearchTask");
    }

    public ParallelSearchTask(List<DbForTask> transientDatabases) : base(MyTask.ParallelSearch)
    {
        // Initialize with appropriate defaults
        SearchParameters = new ParallelSearchParameters()
        {
            TransientDatabases = transientDatabases
        };
        CommonParameters = new(taskDescriptor: "ParallelSearchTask");
    }

    // Rename the TOML section for ParallelSearchParameters to avoid conflicts
    [TomlIgnore]
    public override SearchParameters SearchParameters
    {
        get => ParallelSearchParameters;
        set
        {
            if (value is ParallelSearchParameters msp)
                ParallelSearchParameters = msp;
            else
            {
                // If someone tries to set a base SearchParameters, convert it
                ParallelSearchParameters = new ParallelSearchParameters(value);
            }
        }
    }

    public ParallelSearchParameters ParallelSearchParameters { get; set; } = new();

    #region Properties that are loaded once during Initialization

    [TomlIgnore] public List<DbForTask> PersistentDatabases { get; private set; } = [];

    // Normal Search Stuff
    [TomlIgnore] public MyFileManager MyFileManager = null!;
    [TomlIgnore] public List<Modification> VariableModifications { get; private set; } = [];
    [TomlIgnore] public List<Modification> FixedModifications { get; private set; } = [];
    [TomlIgnore] public List<string> LocalizableModificationTypes { get; private set; } = [];
    [TomlIgnore] public List<IBioPolymer> BaseBioPolymers { get; private set; } = [];
    [TomlIgnore] public Ms2ScanWithSpecificMass[] AllSortedMs2Scans { get; private set; } = [];
    [TomlIgnore] private SpectralMatch[] BaseSearchPsms = null!; // PSMs from base database search
    [TomlIgnore] private PerFileTolerances? _perFileTolerances; // null when every file uses the task's tolerances


    // Optimization caches for FDR alignment and parsimony
    [TomlIgnore] private BaselineScoreLookup<FdrInfo> _baselinePsmFdr = new([]);
    [TomlIgnore] private BaselineScoreLookup<FdrInfo> _baselinePeptideFdr = new([]);
    [TomlIgnore] private ProteinGroupFdrLookup _baselineProteinGroupFdr = new([], FilterType.QValue);
    [TomlIgnore] private List<CachedProteinGroup> CachedBaselineProteinGroups { get; set; } = [];
    [TomlIgnore] private readonly Dictionary<string, List<int>> _baselinePeptideKeyToScanIndexes = new(StringComparer.Ordinal);
    [TomlIgnore] private Dictionary<string, IrtCalibrationModel> _retentionTimeCalibrations = new();

    // Properties to track progress and results across databases
    [TomlIgnore] public int TotalDatabases => ParallelSearchParameters.TransientDatabases.Count;
    [TomlIgnore] public int TotalMs2Scans => AllSortedMs2Scans.Length;
    [TomlIgnore] public int PersistentDatabasePeptideCount { get; set; } = -1; // populated during initialization for use in analysis

    #endregion

    /// <summary> Rough memory for one database searched at once: its copy of the baseline PSM array, its per-scan lock and any PSMs it changes. </summary>
    private const int TransientBytesPerScan = 128;

    protected override MyTaskResults RunSpecific(string outputFolder,
        List<DbForTask> dbFilenameList, List<string> currentRawFileList,
        string taskId, FileSpecificParameters[] fileSettingsList)
    {
        ValidateSettings(fileSettingsList);

        MyTaskResults = new MyTaskResults(this);
        PersistentDatabases = dbFilenameList;

        if (currentRawFileList.All(p => p.EndsWith(".mgf", StringComparison.OrdinalIgnoreCase) || p.EndsWith(".msAlign", StringComparison.OrdinalIgnoreCase)))
        {
            CommonParameters.DoPrecursorDeconvolution = false;
            CommonParameters.UseProvidedPrecursorInfo = true;
        }

        // Initialize unified results manager
        _resultsManager = CreateResultsManager(outputFolder);
        if (_resultsManager.UnreadableSummaryMessage is { } unreadableSummaryMessage)
        {
            Warn(unreadableSummaryMessage);
        }

        // Check cache status early for fast-path optimization
        var allDatabaseNames = ParallelSearchParameters.TransientDatabases
            .Select(db => Path.GetFileNameWithoutExtension(db.FilePath))
            .ToList();
         
        var cacheSummary = _resultsManager!.GetCacheSummary(allDatabaseNames);
        Status(cacheSummary.ToString(), taskId);
        InitializeProgressCounts(cacheSummary);

        // Fast path: If all databases are cached and not overwriting, skip to finalization
        if (cacheSummary.DatabasesNeedingProcessing == 0 && !ParallelSearchParameters.OverwriteTransientSearchOutputs)
        {
            Status("All databases cached, skipping search phase and proceeding to finalization...", taskId);
            goto Finalization;
        }

        // Initialize all necessary data structures including base search
        Initialize(taskId, dbFilenameList, currentRawFileList, fileSettingsList, outputFolder);
        InitializeCompletedDatabaseWriter();

        Status($"Starting search of {TotalDatabases} transient databases...", taskId);

        // The user's MaxThreadsToUsePerFile is one budget shared by every database searched at once; MaxSearchesInParallel caps how many.
        // Each worker gets its own copy of the parameters, so the task's own CommonParameters are never changed.
        int databasesToProcess = Math.Max(1, ParallelSearchParameters.TransientDatabases.Count(db =>
            ParallelSearchParameters.OverwriteTransientSearchOutputs
            || !_resultsManager!.HasCachedResults(Path.GetFileNameWithoutExtension(db.FilePath))));
        long bytesPerDatabase = (long)AllSortedMs2Scans.Length * TransientBytesPerScan;
        FileParallelismPlan plan = FileParallelism.Decide(databasesToProcess, CommonParameters.MaxThreadsToUsePerFile,
            IndexPartitioning.AvailableBytes(), bytesPerDatabase, ParallelSearchParameters.MaxSearchesInParallel);
        ProseCreatedWhileRunning.Append($"transient databases searched in parallel = {plan.FilesInParallel} of {databasesToProcess} " +
            $"(limited by {plan.LimitedBy}), each with {plan.ThreadsPerFile} threads; ");
        CommonParameters transientParameters = CommonParameters.CloneWithNewMaxThreadsToUsePerFile(plan.ThreadsPerFile);

        // Loop through each transient database
        try
        {
            var databases = ParallelSearchParameters.TransientDatabases;
            FileParallelism.ForEachFile(databases.Count, plan.FilesInParallel,
                index => ProcessTransientDatabase(databases[index], outputFolder, taskId, transientParameters));
        }
        finally
        {
            CompleteCompletedDatabaseWriter();
        }

          Finalization:

        // Rewrite the summary so it holds every database, including rows carried over from an earlier run
        Status("Writing summary...", taskId);
        _resultsManager!.PersistAnalysisCache();
        FinishedWritingFile(_resultsManager.SummaryFilePath, [taskId]);

        MyTaskResults.AddTaskSummaryText(
            $"Searched {TotalDatabases} transient databases against {currentRawFileList.Count} spectra files.");

        ReportProgress(new(100, "Parallel search task complete!", [taskId]));
        return MyTaskResults;
    }

    #region Initialization

    /// <summary>
    /// Rejects settings the parallel search cannot honour, before anything is loaded or written.
    /// The transient search engine skips XCorr scoring for speed, and the
    /// summary cache and per-database output folders are keyed by database file name.
    /// </summary>
    private void ValidateSettings(FileSpecificParameters?[]? fileSettingsList)
    {
        var transientDatabases = ParallelSearchParameters.TransientDatabases;
        if (transientDatabases == null || transientDatabases.Count == 0)
        {
            throw new MetaMorpheusException("Parallel search needs at least one transient database. "
                + "Add transient databases to the task before running it.");
        }

        var duplicateNames = transientDatabases
            .GroupBy(db => Path.GetFileNameWithoutExtension(db.FilePath), StringComparer.OrdinalIgnoreCase)
            .Where(g => g.Count() > 1)
            .Select(g => $"{g.Key} ({string.Join(", ", g.Select(db => db.FilePath))})")
            .ToList();
        if (duplicateNames.Count > 0)
        {
            throw new MetaMorpheusException("Parallel search needs every transient database to have a different file name, "
                + "because results are stored by name. These names are used more than once: "
                + string.Join("; ", duplicateNames));
        }

        if (GlobalVariables.AnalyteType == AnalyteType.Oligo)
        {
            throw new MetaMorpheusException("Parallel search does not support RNA searches yet.");
        }

        ValidateDissociationType(CommonParameters.DissociationType, "the task settings");
        if (fileSettingsList != null)
        {
            foreach (var fileSettings in fileSettingsList)
            {
                if (fileSettings?.DissociationType is { } fileDissociationType)
                {
                    ValidateDissociationType(fileDissociationType, "a file-specific settings file");
                }
            }
        }

        if (SearchParameters.SearchType != SearchType.Classic)
        {
            throw new MetaMorpheusException($"Parallel search only runs a classic search, but the search type is set to {SearchParameters.SearchType}. "
                + "Set the search type to Classic and run again.");
        }

        if (SearchParameters.MassDiffAcceptorType.IsMostAbundant()
            || CommonParameters.PrecursorMassMatchMode == PrecursorMassMatchMode.MostAbundant)
        {
            throw new MetaMorpheusException("Parallel search does not support most-abundant precursor matching. "
                + "Choose a monoisotopic precursor mass tolerance and run again.");
        }
    }

    private static void ValidateDissociationType(DissociationType dissociationType, string source)
    {
        if (dissociationType == DissociationType.Autodetect)
        {
            throw new MetaMorpheusException($"Parallel search does not support autodetecting the dissociation type (set in {source}). "
                + "Choose a dissociation type and run again.");
        }

        if (dissociationType == DissociationType.LowCID)
        {
            throw new MetaMorpheusException($"Parallel search does not support Low CID dissociation (set in {source}), "
                + "because it does not score with XCorr. Choose another dissociation type and run again.");
        }
    }

    private void Initialize(string taskId, List<DbForTask> dbFilenameList,
        List<string> currentRawFileList, FileSpecificParameters[] fileSettingsList,
        string outputFolder)
    {
        // Initialize base objects
        MyFileManager = new MyFileManager(SearchParameters.DisposeOfFileWhenDone);

        Status("Loading modifications...", taskId);

        // 1. Load modifications once
        LoadModifications(taskId, out var variableModifications,
            out var fixedModifications, out var localizableModificationTypes);
        VariableModifications = variableModifications;
        FixedModifications = fixedModifications;
        LocalizableModificationTypes = localizableModificationTypes;

        Status("Loading base database(s)...", taskId);

        // 2. Load base database(s) once
        var baseDbLoader = new DatabaseLoadingEngine(CommonParameters,
            FileSpecificParameters, [taskId], dbFilenameList, taskId,
            SearchParameters.DecoyType, SearchParameters.SearchTarget,
            LocalizableModificationTypes);
        BaseBioPolymers = (baseDbLoader.Run() as DatabaseLoadingEngineResults)!.BioPolymers;

        Status($"Loaded {BaseBioPolymers.Count} base proteins", taskId);

        // 3. Load all spectra files once and store in memory
        Status("Loading spectra files...", taskId);
        ConcurrentDictionary<string, Ms2ScanWithSpecificMass[]> loadedSpectraByFile = new();
        int totalMs2Scans = LoadSpectraFiles(currentRawFileList, fileSettingsList, MyFileManager,
            loadedSpectraByFile, taskId);
        AllSortedMs2Scans = loadedSpectraByFile
            .SelectMany(p => p.Value)
            .OrderBy(b => b.PrecursorMass)
            .ToArray();

        // Calibrated per-file tolerances: each scan carries its file's index so the one merged search can honour them
        _perFileTolerances = PerFileTolerances.Create(AllSortedMs2Scans, FileSpecificParameters, CommonParameters, out string? toleranceWarning);
        if (toleranceWarning != null)
            Warn(toleranceWarning);

        // 4. Perform base database search once and store results
        Status("Performing base database search...", taskId);
        BaseSearchPsms = new SpectralMatch[AllSortedMs2Scans.Length];
        PerformSearch(BaseBioPolymers, BaseSearchPsms, [taskId], out _, out int persistentDatabasePeptideCount);
        PersistentDatabasePeptideCount = persistentDatabasePeptideCount;
        Status($"Base search complete. Found {BaseSearchPsms.Count(p => p != null)} PSMs.", taskId);

        // 5. Run baseline FDR once and cache score->q-value mapping for transient alignment
        var baselinePsms = BaseSearchPsms.Where(p => p != null)
            .OrderByDescending(p => p).ToList();

        if (baselinePsms.Count > 0)
        {
            int numNotches = GetMassDiffAcceptor(CommonParameters.PrecursorMassTolerance,
                SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac).NumNotches;
            var baselineFdrEngine = new FdrAnalysisEngine(
                baselinePsms, numNotches, CommonParameters,
                FileSpecificParameters, [taskId], "PSM", false, outputFolder);
            baselineFdrEngine.Run();
        }

        _baselinePsmFdr = BaselineFdr.CreatePsmLookup(baselinePsms);
        _baselinePeptideFdr = BaselineFdr.CreatePeptideLookup(baselinePsms);
        BuildBaselinePeptideToScanIndexLookup();

        // Map Chronologer's % acetonitrile onto each file's minutes once, so every database's RT errors share one scale
        Status("Calibrating predicted retention times...", taskId);
        try
        {
            _retentionTimeCalibrations = RetentionTimeCollector.FitRetentionTimeCalibrations(baselinePsms,
                Math.Min(CommonParameters.QValueThreshold, CommonParameters.PepQValueThreshold),
                CommonParameters.MaxThreadsToUsePerFile, Warn);
        }
        catch (Exception e)
        {
            // e.g. TorchSharp's native libraries are missing; a missing predictor costs the RT columns, not the search
            Warn($"Retention time errors are not reported: Chronologer could not run ({e.GetBaseException().Message})");
        }

        if (SearchParameters.DoParsimony)
        {
            Status("Preparing baseline parsimony cache...", taskId);

            var baselinePsmsForParsimony = FilteredPsms.Filter(baselinePsms,
                commonParams: CommonParameters,
                includeDecoys: true,
                includeContaminants: true,
                includeAmbiguous: false,
                includeHighQValuePsms: true);

            ProteinParsimonyResults baselineParsimonyResults = (ProteinParsimonyResults)new ProteinParsimonyEngine(
                baselinePsmsForParsimony, SearchParameters.ModPeptidesAreDifferent,
                CommonParameters, FileSpecificParameters, [taskId]).Run();

            if (baselineParsimonyResults.ProteinGroups.Count > 0)
            {
                var baselineProteinScoringEngine = new ProteinScoringAndFdrEngine(
                        baselineParsimonyResults.ProteinGroups,
                        baselinePsms,
                        psmsAlreadyFiltered: true, // score every baseline PSM; transient FDR is aligned to it
                        SearchParameters.NoOneHitWonders,
                        SearchParameters.ModPeptidesAreDifferent,
                        true,
                        CommonParameters,
                        FileSpecificParameters,
                        [taskId]);
                var baselineProteinScoringResults = (ProteinScoringAndFdrResults)baselineProteinScoringEngine.Run();

                CachedBaselineProteinGroups = baselineProteinScoringResults.SortedAndScoredProteinGroups
                    .Select(p => new CachedProteinGroup(p))
                    .ToList();

                _baselineProteinGroupFdr = BaselineFdr.CreateProteinGroupLookup(baselineProteinScoringResults
                    .SortedAndScoredProteinGroups, baselineProteinScoringEngine.FilterType);
            }
            else
            {
                CachedBaselineProteinGroups = [];
                _baselineProteinGroupFdr = new([], FilterType.QValue);
            }
        }

        // Write Psms and Peptides from base search
        string psmFile = Path.Combine(outputFolder,
            $"Base{GlobalVariables.AnalyteType.GetSpectralMatchLabel()}s.{GlobalVariables.AnalyteType.GetSpectralMatchExtension()}");
        var psmTask = WritePsmsToTsvAsync(BaseSearchPsms.Where(p => p != null), psmFile, false)
            .ContinueWith(_ => FinishedWritingFile(psmFile, [taskId]));

        var allPeptides = BaseSearchPsms.Where(p => p != null).CollapseToPeptides(true).OrderByDescending(p => p).ToList();
        string peptideFile = Path.Combine(outputFolder,
            $"Base{GlobalVariables.AnalyteType}s.{GlobalVariables.AnalyteType.GetSpectralMatchExtension()}");
        var peptideTask = WritePsmsToTsvAsync(allPeptides, peptideFile, true)
            .ContinueWith(_ => FinishedWritingFile(peptideFile, [taskId]));

        Task proteinTask = Task.CompletedTask;
        if (SearchParameters.DoParsimony && CachedBaselineProteinGroups.Any())
        {
            var groups = CachedBaselineProteinGroups.Select(p => {
                var group = p.GetProteinGroup();
                group.GetIdentifiedPeptidesOutput(SearchParameters.SilacLabels);
                return group;                
                });
            string proteinFile = Path.Combine(outputFolder,
                $"Base{GlobalVariables.AnalyteType.GetBioPolymerLabel()}Groups.tsv");
            proteinTask = WriteProteinGroupsToTsvAsync(groups.ToList(), proteinFile)
                .ContinueWith(_ => FinishedWritingFile(proteinFile, [taskId]));
        }

        // Write prose for base settings
        ProseCreatedWhileRunning.Append(
            $"Base database contained {BaseBioPolymers.Count(p => !p.IsDecoy)} non-decoy protein entries. ");
        ProseCreatedWhileRunning.Append(
            $"Searching {ParallelSearchParameters.TransientDatabases.Count} transient databases against {currentRawFileList.Count} spectra files. ");

        Task.WhenAll([psmTask, peptideTask, proteinTask]).Wait();
    }

    private int LoadSpectraFiles(List<string> currentRawFileList, FileSpecificParameters[] fileSettingsList,
    MyFileManager myFileManager, ConcurrentDictionary<string, Ms2ScanWithSpecificMass[]> loadedSpectraByFile,
    string taskId)
    {
         int totalMs2Scans = 0;
         int specLoadingProgress = 0;
         var specLoadingNestedIds = new List<string> { taskId, "Spectra Loading" };
         Status("Loading spectra files...", specLoadingNestedIds);

        // The user's thread budget is shared by the files loaded at once; copies are made here, on one thread.
        FileParallelismPlan plan = FileParallelism.Decide(currentRawFileList.Count, CommonParameters.MaxThreadsToUsePerFile,
            0, 0, ParallelSearchParameters.MaxSearchesInParallel);
        var paramsPerFile = currentRawFileList
            .Select((rawFile, i) => SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[i], rawFile)
                .CloneWithNewMaxThreadsToUsePerFile(plan.ThreadsPerFile))
            .ToArray();

        FileParallelism.ForEachFile(currentRawFileList.Count, plan.FilesInParallel,
             fileIndex =>
             {
                 string rawFile = currentRawFileList[fileIndex];
                 var fileParams = paramsPerFile[fileIndex];
                 var msDataFile = myFileManager.LoadFile(rawFile, fileParams);
                 var ms2Scans = GetMs2Scans(msDataFile, rawFile, fileParams)
                     .OrderBy(b => b.PrecursorMass).ToArray();
                 loadedSpectraByFile.AddOrUpdate(rawFile, ms2Scans, (key, oldValue) => ms2Scans);
                 Interlocked.Add(ref totalMs2Scans, ms2Scans.Length);
                 myFileManager.DoneWithFile(rawFile);

                 lock (_progressLock)
                 {
                     ReportProgress(new ProgressEventArgs(
                         (int)(Interlocked.Increment(ref specLoadingProgress) / (double)currentRawFileList.Count * 100),
                         $"Loaded {Path.GetFileName(rawFile)}",
                         specLoadingNestedIds));
                 }
             });

         ReportProgress(new ProgressEventArgs(100, $"Finished Loading spectra files.", specLoadingNestedIds));
         Status($"Finished Loading {currentRawFileList.Count} spectra files.", taskId);

         return totalMs2Scans;
    }

    /// <summary>
    /// Creates the results manager with the configured metric collectors
    /// </summary>
    private TransientDatabaseResultsManager CreateResultsManager(string outputFolder)
    {
        // Define cache paths
        string analysisCachePath = Path.Combine(outputFolder, TransientDatabaseResultsManager.SummaryResultsFileName);

        // Initialize collectors based on user preferences
        var collectors = new List<IMetricCollector>
        {
            new BasicMetricCollector(),
            new PsmPeptideSearchCollector("Homo sapiens"),
            new FragmentIonCollector(),
            new RetentionTimeCollector(),
        };

        if (ParallelSearchParameters.DoParsimony)
        {
            collectors.Add(new ProteinGroupCollector("Homo sapiens"));
        }

        return new TransientDatabaseResultsManager(new MetricAggregator(collectors), analysisCachePath);
    }

    #endregion

    #region Search

    private void ProcessTransientDatabase(DbForTask transientDb, string outputFolder, string taskId, CommonParameters commonParameters)
    {
        if (GlobalVariables.StopLoops)
            return;

         string dbName = Path.GetFileNameWithoutExtension(transientDb.FilePath);
         string dbOutputFolder = Path.Combine(outputFolder, dbName);
         List<string> nestedIds = [taskId, dbName];

        Status($"Processing {dbName}...", nestedIds);

         // Check if we should skip or overwrite this database
         bool shouldProcess = !_resultsManager!.HasCachedResults(dbName) || 
                            ParallelSearchParameters.OverwriteTransientSearchOutputs;

        if (!shouldProcess)
        {
            ReportProgress(new(100, $"Skipping {dbName} - results already exist in cache", nestedIds));
            UpdateProgress(TotalDatabases, taskId);
            return;
        }


         // Handle overwrite scenario
        if (ParallelSearchParameters.OverwriteTransientSearchOutputs && _resultsManager.HasCachedResults(dbName))
        {
            Status($"Overwriting existing results for {dbName}...", nestedIds);
            if (Directory.Exists(dbOutputFolder))
            {
                Directory.Delete(dbOutputFolder, true);
            }
        }

         if (!Directory.Exists(dbOutputFolder))
             Directory.CreateDirectory(dbOutputFolder);

        Status($"Loading transient database {dbName}...", nestedIds);

         // Load transient database
         var transientProteins = LoadTransientDatabase(transientDb, nestedIds, taskId, commonParameters);

         if (GlobalVariables.StopLoops)
         {
             return;
         }

         // Create HashSet of transient bioPolymer accessions for later filtering
         var transientProteinAccessions = new HashSet<string>(
             transientProteins.Select(p => p.Accession));

        Status($"Searching {dbName} ({transientProteins.Count} transient proteins)...", nestedIds);

         // Reuse baseline PSMs with copy-on-write in peptide/proteoform mode.
         SpectralMatch[] psmArray = BaseSearchPsms.ToArray();
         PerformSearch(transientProteins, psmArray, nestedIds, out HashSet<int> updatedPsmIndexes, out int transientPeptideCount, useCopyOnWrite: true, commonParameters: commonParameters);

        Status($"Performing post-search analysis for {dbName}...", nestedIds);

         int totalProteins = BaseBioPolymers.Count + transientProteins.Count;
         
         // Process database through unified manager (handles analysis + caching)
         var (analysisContext, dbResults) = PerformPostSearchAnalysis(
              psmArray,
              dbOutputFolder, 
              nestedIds,
             dbName, 
             totalProteins, 
             transientProteinAccessions, 
             transientPeptideCount, 
             transientProteins, 
              transientDb,
              updatedPsmIndexes,
              commonParameters
          );

         _completedDatabaseWriteChannel!.Writer.WriteAsync((analysisContext, dbResults)).AsTask().GetAwaiter().GetResult();

        // Cleanup transient proteins to free memory
        transientProteins.Clear();

        ReportProgress(new(100, $"Finished analysis for {dbName}; queued output writing", nestedIds));
    }

    /// <summary>
    /// Populates and returns the spectral match array using classic search engine
    /// </summary>
     private void PerformSearch(List<IBioPolymer> proteinsToSearch, SpectralMatch[] spectralMatchArray, List<string> nestedIds, out HashSet<int> updatedPsmIndexes, out int peptidesSearched, bool useCopyOnWrite = false, CommonParameters? commonParameters = null)
     {
         commonParameters ??= CommonParameters;
         // with per-file tolerances the window is searched at its widest and narrowed per scan inside the engine
         var massDiffAcceptor = GetMassDiffAcceptor(
             _perFileTolerances?.WidestPrecursor ?? CommonParameters.PrecursorMassTolerance,
             SearchParameters.MassDiffAcceptorType,
             SearchParameters.CustomMdac);

         // Run the classic search engine
          var searchEngine = new TransientClassicSearchEngine(
              spectralMatchArray, AllSortedMs2Scans, VariableModifications,
              FixedModifications, proteinsToSearch, massDiffAcceptor, commonParameters,
              FileSpecificParameters, nestedIds, copyOnWriteEnabled: useCopyOnWrite,
              perFileTolerances: _perFileTolerances);

         var results = searchEngine.Run();
         updatedPsmIndexes = (results as TransientSearchEngineResults)!.UpdatedSpectralMatchIndexes;
         peptidesSearched = (results as TransientSearchEngineResults)!.PeptidesSearched;
         ReportProgress(new(100, "Finished Classic Search...", nestedIds));
     }

    private (TransientDatabaseContext Context, TransientDatabaseMetrics Metrics) PerformPostSearchAnalysis(
        SpectralMatch[] allPsms,
        string outputFolder,
        List<string> nestedIds,
        string dbName,
        int totalProteins,
        HashSet<string> transientProteinAccessions,
        int transientPeptideCount,
        List<IBioPolymer> transientProteins,
        DbForTask transientDatabase,
        HashSet<int> updatedPsmIndexes,
        CommonParameters commonParameters)
    {
        bool writeAllResults = !ParallelSearchParameters.WriteTransientResultsOnly;

        List<SpectralMatch> psmList = allPsms.Where(p => p != null)
            .OrderByDescending(p => p).ToList();

        #region FDR and Parsiomony

        Status($"Performing FDR Analysis for {dbName}...", nestedIds);

        // Disambiguate - modify PSMs in place - Only used for notch disambiguation and parallel search does not use notches. 
        //var disambiguationEngine = new DisambiguationEngine(
        //    psmList, CommonParameters, FileSpecificParameters, nestedIds);
        //await disambiguationEngine.RunAsync();
        //DebugStatus("Disambiguation complete", nestedIds, dbName, postAnalysisStopwatch);

        // Apply FDR Alignment to PSMs
        // Only scans the transient search changed can hold a transient protein, so only those are checked and sorted.
        // Taking them in scan order first keeps ties in the same order as the full sorted list.
        var transientPsms = FilterToTransientDatabaseOnly(
                updatedPsmIndexes.Order().Select(i => allPsms[i]).Where(p => p != null),
                transientProteinAccessions)
            .OrderByDescending(p => p)
            .ToList();
        BaselineFdr.AlignPsms(_baselinePsmFdr, transientPsms);

        // Collapse to peptides and apply FDR alignment to peptides
        var allPeptides = psmList.CollapseToPeptides(true).ToList();
        var transientPeptides = FilterToTransientDatabaseOnly(allPeptides, transientProteinAccessions).ToList();
        BaselineFdr.AlignPeptides(_baselinePeptideFdr, transientPeptides);

        List<ProteinGroup>? proteinGroups = null;
        if (SearchParameters.DoParsimony && transientPsms.Count > 0)
        {
            Status($"Performing parsimony for {dbName}...", nestedIds);

            var transientParsimonyEngine = new TransientProteinParsimonyEngine(
                allPsms,
                BaseSearchPsms,
                updatedPsmIndexes,
                transientProteins,
                _baselinePeptideKeyToScanIndexes,
                SearchParameters.ModPeptidesAreDifferent,
                commonParameters,
                FileSpecificParameters,
                nestedIds);

            ProteinParsimonyResults proteinAnalysisResults =
                (ProteinParsimonyResults)transientParsimonyEngine.Run();
            List<ProteinGroup> transientParsimonyGroups = proteinAnalysisResults.ProteinGroups;
            List<SpectralMatch> transientPsmsForParsimony = transientParsimonyEngine.NeighborhoodPsms;

            //var baselineRuntimeGroups = CachedBaselineProteinGroups
            //    .Select(group => group.CreateRuntimeCopy())
            //    .ToList();

            if (CachedBaselineProteinGroups.Count > 0 || transientParsimonyGroups.Count > 0)
            {
                // Writing is a mutation process for Protein Groups, so if we are writing ALL (not just transient) we need to make a copy.
                // If we are only writing transient results, we can get away with not copying the baseline groups.
                List<ProteinGroup> baselineProteinGroups =  CachedBaselineProteinGroups
                    .Select(pg => writeAllResults ? pg.CreateRuntimeCopy() : pg.GetProteinGroup())
                    .ToList();

                ProteinScoringAndFdrResults proteinScoringAndFdrResults =
                        (ProteinScoringAndFdrResults)new TransientProteinScoringAndFdrEngine(
                            baselineProteinGroups,
                            transientParsimonyGroups,
                            transientPsmsForParsimony,
                            _baselineProteinGroupFdr,
                            SearchParameters.NoOneHitWonders, SearchParameters.ModPeptidesAreDifferent,
                            true, commonParameters, FileSpecificParameters, nestedIds).Run();

                proteinGroups = proteinScoringAndFdrResults.SortedAndScoredProteinGroups;
            }
            else
            {
                proteinGroups = [];
            }
        }

        #endregion

        List<ProteinGroup>? transientProteinGroups = null;
        if (proteinGroups is not null)
        {
            // Count bioPolymer groups that contain at least one transient database bioPolymer
            transientProteinGroups =
                FilterProteinGroupsToTransientDatabaseOnly(proteinGroups, transientProteinAccessions)
                .Select(CachedProteinGroup.DeepCopy)
                .ToList();
        }

        #region Result Caching and Analysis

        // Create analysis context
        var analysisContext = new TransientDatabaseContext
        {
            DatabaseName = dbName,
            TransientDatabase = transientDatabase,
            TransientProteins = transientProteins,
            TransientProteinAccessions = transientProteinAccessions,
            AllPsms = psmList,
            TransientPsms = transientPsms,
            AllPeptides = allPeptides,
            TransientPeptides = transientPeptides,
            ProteinGroups = proteinGroups,
            TransientProteinGroups = transientProteinGroups,
            RetentionTimeCalibrations = _retentionTimeCalibrations,
            CommonParameters = commonParameters,
            TotalProteins = totalProteins,
            TransientPeptideCount = transientPeptideCount,
            OutputFolder = outputFolder,
            NestedIds = nestedIds
        };

        // Process through unified manager (caches analysis results)
        Status($"Running analysis for {dbName}...", nestedIds);
        var aggregatedResult = _resultsManager!.ProcessDatabase(
            analysisContext,
            forceRecompute: ParallelSearchParameters.OverwriteTransientSearchOutputs
        );

        #endregion

        return (analysisContext, aggregatedResult);
    }

    #endregion

    #region Result Writing

    private void InitializeCompletedDatabaseWriter()
    {
        _completedDatabaseWriteChannel = Channel.CreateBounded<(TransientDatabaseContext Context, TransientDatabaseMetrics Metrics)>(
            new BoundedChannelOptions(2)
            {
                SingleReader = true,
                SingleWriter = false,
                FullMode = BoundedChannelFullMode.Wait
            });
        _completedDatabaseWriterTask = Task.Run(RunCompletedDatabaseWriterLoopAsync);
    }

    private void CompleteCompletedDatabaseWriter()
    {
        if (_completedDatabaseWriteChannel == null)
            return;

        _completedDatabaseWriteChannel.Writer.TryComplete();
        _completedDatabaseWriterTask?.GetAwaiter().GetResult();
        _completedDatabaseWriteChannel = null;
        _completedDatabaseWriterTask = null;
    }

    private async Task RunCompletedDatabaseWriterLoopAsync()
    {
        var channel = _completedDatabaseWriteChannel ?? throw new InvalidOperationException("Completed database writer channel is not initialized.");

        try
        {
            await foreach (var (context, metrics) in channel.Reader.ReadAllAsync())
            {
                await WriteCompletedDatabaseOutputsAsync(context, metrics);
            }
        }
        catch (Exception ex)
        {
            channel.Writer.TryComplete(ex);
            throw;
        }
    }

    private async Task WriteCompletedDatabaseOutputsAsync(TransientDatabaseContext context, TransientDatabaseMetrics metrics)
    {
        string dbName = context.DatabaseName;
        string outputFolder = context.OutputFolder;
        List<string> nestedIds = context.NestedIds;
        bool writeAllResults = !ParallelSearchParameters.WriteTransientResultsOnly;

        Status($"Writing results for {dbName}...", nestedIds);

        string transientPsmFile = Path.Combine(outputFolder,
            $"{dbName}_All{GlobalVariables.AnalyteType.GetSpectralMatchLabel()}s.{GlobalVariables.AnalyteType.GetSpectralMatchExtension()}");
        await WritePsmsToTsvAsync(context.TransientPsms, transientPsmFile, false);
        FinishedWritingFile(transientPsmFile, nestedIds);

        if (writeAllResults)
        {
            string psmFile = Path.Combine(outputFolder,
                $"All{GlobalVariables.AnalyteType.GetSpectralMatchLabel()}s.{GlobalVariables.AnalyteType.GetSpectralMatchExtension()}");
            await WritePsmsToTsvAsync(context.AllPsms, psmFile, false);
            FinishedWritingFile(psmFile, nestedIds);
        }

        if (context.TransientPeptides is { Count: > 0 })
        {
            string transientPeptideFile = Path.Combine(outputFolder,
            $"{dbName}_All{GlobalVariables.AnalyteType}s.{GlobalVariables.AnalyteType.GetSpectralMatchExtension()}");
            await WritePsmsToTsvAsync(context.TransientPeptides, transientPeptideFile, true);
            FinishedWritingFile(transientPeptideFile, nestedIds);
        }

        if (writeAllResults)
        {
            string peptideFile = Path.Combine(outputFolder,
                $"All{GlobalVariables.AnalyteType}s.{GlobalVariables.AnalyteType.GetSpectralMatchExtension()}");
            await WritePsmsToTsvAsync(context.AllPeptides, peptideFile, true);
            FinishedWritingFile(peptideFile, nestedIds);
        }

        if (context.TransientProteinGroups is { Count: > 0 })
        {
            context.TransientProteinGroups.ForEach(x => x.GetIdentifiedPeptidesOutput(SearchParameters.SilacLabels));
            string transientProteinFile = Path.Combine(outputFolder,
                $"{dbName}_All{GlobalVariables.AnalyteType.GetBioPolymerLabel()}Groups.tsv");
            await WriteProteinGroupsToTsvAsync(context.TransientProteinGroups, transientProteinFile);
            FinishedWritingFile(transientProteinFile, nestedIds);
        }

        if (writeAllResults && context.ProteinGroups is { Count: > 0 })
        {
            context.ProteinGroups.ForEach(x => x.GetIdentifiedPeptidesOutput(SearchParameters.SilacLabels));
            string proteinFile = Path.Combine(outputFolder,
                $"All{GlobalVariables.AnalyteType.GetBioPolymerLabel()}Groups.tsv");
            await WriteProteinGroupsToTsvAsync(context.ProteinGroups, proteinFile);
            FinishedWritingFile(proteinFile, nestedIds);
        }

        if (SearchParameters.WriteSpectralLibrary)
        {
            string spectralLibraryPath = Path.Combine(outputFolder, $"AllPeptidesAnd_{dbName}_SpectralLibrary.msp");
            await WriteSpectralLibraryAsync(context.AllPsms, spectralLibraryPath);
            FinishedWritingFile(spectralLibraryPath, nestedIds);
        }

        if (ParallelSearchParameters.WriteTransientSpectralLibrary)
        {
            string spectralLibraryPath = Path.Combine(outputFolder, $"{dbName}_SpectralLibrary.msp");
            await WriteSpectralLibraryAsync(context.TransientPsms, spectralLibraryPath);
            FinishedWritingFile(spectralLibraryPath, nestedIds);
        }

        await WriteIndividualDatabaseResultsTextAsync(metrics, outputFolder, nestedIds);

        if (ParallelSearchParameters.CompressTransientSearchOutputs)
        {
            Status($"Compressing output for {dbName}...", nestedIds);
            CompressTransientDatabaseOutput(outputFolder);
        }

        _resultsManager!.AppendCheckpoint(metrics);
        MarkDatabaseCompleted();
        UpdateProgress(TotalDatabases, nestedIds[0]);
        ReportProgress(new(100, $"Finished {dbName}", nestedIds));
    }

    private async Task WriteIndividualDatabaseResultsTextAsync(TransientDatabaseMetrics results, string outputFolder, List<string> nestedIds)
    {
        var resultsPath = Path.Combine(outputFolder, "results.txt");
        await results.WriteToTextFileAsync(resultsPath, CommonParameters.QValueThreshold, SearchParameters.DoParsimony);
        FinishedWritingFile(resultsPath, nestedIds);
    }

    private async Task WriteProteinGroupsToTsvAsync(List<ProteinGroup> proteinGroups, string filePath)
    {
        if (proteinGroups != null && proteinGroups.Any())
        {
            double qValueThreshold = Math.Min(CommonParameters.QValueThreshold, CommonParameters.PepQValueThreshold);
            using (StreamWriter output = new StreamWriter(filePath))
            {
                await output.WriteLineAsync(proteinGroups.First().GetTabSeparatedHeader());
                for (int i = 0; i < proteinGroups.Count; i++)
                {
                    if (!SearchParameters.WriteDecoys && proteinGroups[i].IsDecoy ||
                        !SearchParameters.WriteContaminants && proteinGroups[i].IsContaminant ||
                        !SearchParameters.WriteHighQValuePsms && proteinGroups[i].QValue > qValueThreshold)
                    {
                        continue;
                    }
                    else
                    {
                        await output.WriteLineAsync(proteinGroups[i].ToString());
                    }
                }
            }
        }
    }

    /// <summary>
    /// Writes PSMs or peptides as a normal search does: the WriteDecoys, WriteContaminants and
    /// WriteHighQValuePsms filter from <see cref="PostSearchAnalysisTask"/>, sorted best first, then
    /// <see cref="MetaMorpheusTask.WritePsmsToTsv"/> for the columns (1/K0, collision energy and the
    /// most-abundant mass error when the data has them). Runs off the caller's thread so the baseline's
    /// files can be written at the same time.
    /// </summary>
    private Task WritePsmsToTsvAsync(IEnumerable<SpectralMatch> psms, string filePath, bool writePeptideLevelResults)
    {
        return Task.Run(() =>
        {
            var psmsToWrite = FilteredPsms.Filter(psms,
                    CommonParameters,
                    includeDecoys: SearchParameters.WriteDecoys,
                    includeContaminants: SearchParameters.WriteContaminants,
                    includeAmbiguous: true,
                    includeHighQValuePsms: SearchParameters.WriteHighQValuePsms,
                    filterAtPeptideLevel: writePeptideLevelResults)
                .OrderByDescending(p => p)
                .ToList();

            WritePsmsToTsv(psmsToWrite, filePath, SearchParameters.ModsToWriteSelection, writePeptideLevelResults);
        });
    }

    private async Task WriteSpectralLibraryAsync(IEnumerable<SpectralMatch> psms, string outFilePath)
    {
        try
        {
            var peptidesForSpectralLibrary = FilteredPsms.Filter(psms,
                CommonParameters,
                includeDecoys: false,
                includeContaminants: false,
                includeAmbiguous: false,
                includeHighQValuePsms: false);

            //group psms by peptide and charge, the psms having same sequence and same charge will be in the same group
            IEnumerable<LibrarySpectrum> spectraLibrary = peptidesForSpectralLibrary.GroupBy(p => (p.FullSequence, p.ScanPrecursorCharge))
                .Select(p => p.MaxBy(q => q.Score))
                .Where(p => p != null)
                .Select(p => new LibrarySpectrum(
                    p!.FullSequence,
                    p.ScanPrecursorMonoisotopicPeakMz,
                    p.ScanPrecursorCharge,
                    p.MatchedFragmentIons,
                    p.ScanRetentionTime));

            await using StreamWriter output = new StreamWriter(outFilePath);
            foreach (var x in spectraLibrary)
            {
                await output.WriteLineAsync(x.ToString());
            }
        }
        catch (Exception e)
        {
            EngineCrashed("SpectralLibraryGeneration", e);
        }
    }

    /// <summary>
    /// Compresses the transient database output folder
    /// </summary>
    private void CompressTransientDatabaseOutput(string outputFolder)
    {
        var directoryInfo = new DirectoryInfo(outputFolder);
        foreach (FileInfo fileToCompress in directoryInfo.GetFiles())
        {
            bool compressed = false;
            int maxRetries = 10;
            int retryCount = 0;
            while (!compressed && retryCount < maxRetries)
            {
                try
                {
                    MyFileManager.CompressFile(fileToCompress);
                    compressed = true;
                }
                catch (IOException ex) when (ex is IOException && (ex.HResult & 0xFFFF) == 32) // ERROR_SHARING_VIOLATION
                {
                    // File is being used by another process, wait and retry
                    Thread.Sleep(1000);
                    retryCount++;
                }
                catch (Exception ex)
                {
                    Warn($"Failed to compress file {fileToCompress.FullName}: {ex.Message}");
                    break;
                }
            }
            if (!compressed)
            {
                Warn($"Could not compress file {fileToCompress.FullName} after {maxRetries} retries.");
            }
        }
    }

    #endregion

    #region Transient Protein Handling

    private List<IBioPolymer> LoadTransientDatabase(DbForTask transientDbPath,
        List<string> nestedIds, string taskId, CommonParameters commonParameters)
    {
        var transientDbList = new List<DbForTask> { transientDbPath };
        var transientDbLoader = new DatabaseLoadingEngine(commonParameters,
            FileSpecificParameters, nestedIds, transientDbList, taskId,
            SearchParameters.DecoyType, SearchParameters.SearchTarget,
            LocalizableModificationTypes);
        var transientProteins = (transientDbLoader.Run() as DatabaseLoadingEngineResults)!.BioPolymers;
        return transientProteins;
    }

    /// <summary>
    /// Filters spectral matches to only include those that match exclusively to transient database proteins
    /// </summary>
    private IEnumerable<SpectralMatch> FilterToTransientDatabaseOnly(IEnumerable<SpectralMatch> spectralMatches,
        HashSet<string> transientProteinAccessions, bool allowAmbiguous = true)
    {
        foreach (var psm in spectralMatches)
        {
            bool hasTransientProtein = false;

            if (allowAmbiguous)
            {
                foreach (var match in psm.BestMatchingBioPolymersWithSetMods)
                {
                    if (transientProteinAccessions.Contains(match.SpecificBioPolymer.Parent.Accession))
                    {
                        hasTransientProtein = true;
                        break;
                    }
                }
            }
            else if (psm.BestMatchingBioPolymersWithSetMods.All(p => transientProteinAccessions.Contains(p.SpecificBioPolymer.Parent.Accession)))
                hasTransientProtein = true;

            if (hasTransientProtein)
            {
                yield return psm;
            }
        }
    }

    /// <summary>
    /// Filters bioPolymer groups to only include those where all proteins are from the transient database
    /// </summary>
    private IEnumerable<ProteinGroup> FilterProteinGroupsToTransientDatabaseOnly(List<ProteinGroup> proteinGroups,
        HashSet<string> transientProteinAccessions)
    {
        return proteinGroups
            .Where(pg => pg.Proteins.Any(p => transientProteinAccessions.Contains(p.Accession)));
    }

    private void BuildBaselinePeptideToScanIndexLookup()
    {
        _baselinePeptideKeyToScanIndexes.Clear();

        for (int scanIndex = 0; scanIndex < BaseSearchPsms.Length; scanIndex++)
        {
            var psm = BaseSearchPsms[scanIndex];
            if (psm is null)
                continue;

            HashSet<string> keysSeenForScan = new(StringComparer.Ordinal);
            foreach (var hypothesis in psm.BestMatchingBioPolymersWithSetMods)
            {
                string peptideKey = TransientProteinParsimonyEngine.GetParsimonyPeptideKey(
                    hypothesis.SpecificBioPolymer,
                    SearchParameters.ModPeptidesAreDifferent);
                if (!keysSeenForScan.Add(peptideKey))
                    continue;

                if (!_baselinePeptideKeyToScanIndexes.TryGetValue(peptideKey, out var scanIndexes))
                {
                    scanIndexes = [];
                    _baselinePeptideKeyToScanIndexes.Add(peptideKey, scanIndexes);
                }

                scanIndexes.Add(scanIndex);
            }
        }
    }

    #endregion

    #region Progress
    private void InitializeProgressCounts(CacheSummary cacheSummary)
    {
        lock (_countLock)
        {
            _cachedAtStart = cacheSummary.CachedDatabases;
            _completedThisRun = 0;
            _initialFinishedCount = ParallelSearchParameters.OverwriteTransientSearchOutputs
                ? 0
                : cacheSummary.CachedDatabases;
        }
    }

    private void MarkDatabaseCompleted()
    {
        lock (_countLock)
        {
            _completedThisRun++;
        }
    }

    private (int Finished, int Todo, int Cached) GetProgressCounts()
    {
        lock (_countLock)
        {
            int finished = _initialFinishedCount + _completedThisRun;
            finished = Math.Min(finished, TotalDatabases);
            int todo = Math.Max(0, TotalDatabases - finished);
            return (finished, todo, _cachedAtStart);
        }
    }

    private void UpdateProgress(int totalDatabases, string taskId)
    {
        lock (_progressLock)
        {
            var counts = GetProgressCounts();
            ReportProgress(new ProgressEventArgs(
                totalDatabases == 0 ? 100 : (int)(counts.Finished / (double)totalDatabases * 100),
                $"Completed {counts.Finished}/{totalDatabases} databases",
                new List<string> { taskId }));
        }
    }

    #endregion
}
