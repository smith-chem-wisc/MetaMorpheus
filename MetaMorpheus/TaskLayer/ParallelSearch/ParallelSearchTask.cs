#nullable enable
using EngineLayer;
using EngineLayer.DatabaseLoading;
using EngineLayer.FdrAnalysis;
using EngineLayer.ParallelSearch;
using EngineLayer.ParallelSearch.FdrAlignment;
using EngineLayer.SpectrumMatch;
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
                ParallelSearchParameters = new ParallelSearchParameters(SearchParameters);
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


    // Optimization caches for FDR alignment and parsimony
    [TomlIgnore] private BaselineScoreLookup<FdrInfo> _baselinePsmFdr = new([]);
    [TomlIgnore] private BaselineScoreLookup<FdrInfo> _baselinePeptideFdr = new([]);
    [TomlIgnore] private BaselineScoreLookup<ProteinGroupFdr> _baselineProteinGroupFdr = new([]);
    [TomlIgnore] private List<CachedProteinGroup> CachedBaselineProteinGroups { get; set; } = [];
    [TomlIgnore] private readonly Dictionary<string, List<int>> _baselinePeptideKeyToScanIndexes = new(StringComparer.Ordinal);

    // Properties to track progress and results across databases
    [TomlIgnore] public int TotalDatabases => ParallelSearchParameters.TransientDatabases.Count;
    [TomlIgnore] public int TotalMs2Scans => AllSortedMs2Scans.Length;
    [TomlIgnore] public int PersistentDatabasePeptideCount { get; set; } = -1; // populated during initialization for use in analysis

    #endregion

    protected override MyTaskResults RunSpecific(string outputFolder,
        List<DbForTask> dbFilenameList, List<string> currentRawFileList,
        string taskId, FileSpecificParameters[] fileSettingsList)
    {
        MyTaskResults = new MyTaskResults(this);
        PersistentDatabases = dbFilenameList;

        if (currentRawFileList.All(p => p.EndsWith(".mgf", StringComparison.OrdinalIgnoreCase) || p.EndsWith(".msAlign", StringComparison.OrdinalIgnoreCase)))
        {
            CommonParameters.DoPrecursorDeconvolution = false;
            CommonParameters.UseProvidedPrecursorInfo = true;
        }

        // Initialize unified results manager
        _resultsManager = CreateResultsManager(outputFolder);
        
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

        // Determine optimal thread allocation
        int totalAvailableThreads = Environment.ProcessorCount;
         int databaseParallelism = Math.Min(ParallelSearchParameters.MaxSearchesInParallel,
             ParallelSearchParameters.TransientDatabases.Count);
         int threadsPerDatabase = Math.Max(1, totalAvailableThreads / databaseParallelism);
         CommonParameters.MaxThreadsToUsePerFile = threadsPerDatabase;

        // Loop through each transient database
         try
         {
             Parallel.ForEach(ParallelSearchParameters.TransientDatabases,
                 new ParallelOptions { MaxDegreeOfParallelism = databaseParallelism },
                 transientDbPath =>
                 {
                     ProcessTransientDatabase(transientDbPath, outputFolder, taskId);
                 });
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
                ProteinScoringAndFdrResults baselineProteinScoringResults =
                    (ProteinScoringAndFdrResults)new ProteinScoringAndFdrEngine(
                        baselineParsimonyResults.ProteinGroups,
                        baselinePsms,
                        psmsAlreadyFiltered: true, // score every baseline PSM; transient FDR is aligned to it
                        SearchParameters.NoOneHitWonders,
                        SearchParameters.ModPeptidesAreDifferent,
                        true,
                        CommonParameters,
                        FileSpecificParameters,
                        [taskId]).Run();

                CachedBaselineProteinGroups = baselineProteinScoringResults.SortedAndScoredProteinGroups
                    .Select(p => new CachedProteinGroup(p))
                    .ToList();

                _baselineProteinGroupFdr = BaselineFdr.CreateProteinGroupLookup(baselineProteinScoringResults
                    .SortedAndScoredProteinGroups);
            }
            else
            {
                CachedBaselineProteinGroups = [];
                _baselineProteinGroupFdr = new([]);
            }
        }

        // Write Psms and Peptides from base search
        string psmFile = Path.Combine(outputFolder,
            $"Base{GlobalVariables.AnalyteType.GetSpectralMatchLabel()}s.{GlobalVariables.AnalyteType.GetSpectralMatchExtension()}");
        var psmTask = WritePsmsToTsvAsync(BaseSearchPsms.Where(p => p != null).OrderByDescending(p => p), psmFile, SearchParameters.ModsToWriteSelection, false)
            .ContinueWith(_ => FinishedWritingFile(psmFile, [taskId]));

        var allPeptides = BaseSearchPsms.Where(p => p != null).CollapseToPeptides(true).OrderByDescending(p => p).ToList();
        string peptideFile = Path.Combine(outputFolder,
            $"Base{GlobalVariables.AnalyteType}s.{GlobalVariables.AnalyteType.GetSpectralMatchExtension()}");
        var peptideTask = WritePsmsToTsvAsync(allPeptides, peptideFile, SearchParameters.ModsToWriteSelection, true)
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

        Parallel.ForEach(currentRawFileList,
            new ParallelOptions { MaxDegreeOfParallelism = ParallelSearchParameters.MaxSearchesInParallel },
             rawFile =>
             {
                 var fileParams = SetAllFileSpecificCommonParams(CommonParameters,
                     fileSettingsList[currentRawFileList.IndexOf(rawFile)], rawFile);
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

    private void ProcessTransientDatabase(DbForTask transientDb, string outputFolder, string taskId)
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
         var transientProteins = LoadTransientDatabase(transientDb, nestedIds, taskId);

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
         PerformSearch(transientProteins, psmArray, nestedIds, out HashSet<int> updatedPsmIndexes, out int transientPeptideCount, useCopyOnWrite: true);

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
              updatedPsmIndexes
          ).GetAwaiter().GetResult();

         _completedDatabaseWriteChannel!.Writer.WriteAsync((analysisContext, dbResults)).AsTask().GetAwaiter().GetResult();

        // Cleanup transient proteins to free memory
        transientProteins.Clear();

        ReportProgress(new(100, $"Finished analysis for {dbName}; queued output writing", nestedIds));
    }

    /// <summary>
    /// Populates and returns the spectral match array using classic search engine
    /// </summary>
     private void PerformSearch(List<IBioPolymer> proteinsToSearch, SpectralMatch[] spectralMatchArray, List<string> nestedIds, out HashSet<int> updatedPsmIndexes, out int peptidesSearched, bool useCopyOnWrite = false)
     {
         var massDiffAcceptor = GetMassDiffAcceptor(
             CommonParameters.PrecursorMassTolerance,
             SearchParameters.MassDiffAcceptorType,
             SearchParameters.CustomMdac);

         // Run the classic search engine
          var searchEngine = new TransientClassicSearchEngine(
              spectralMatchArray, AllSortedMs2Scans, VariableModifications,
              FixedModifications, proteinsToSearch, massDiffAcceptor, CommonParameters,
              FileSpecificParameters, nestedIds, copyOnWriteEnabled: useCopyOnWrite);

         var results = searchEngine.Run();
         updatedPsmIndexes = (results as TransientSearchEngineResults)!.UpdatedSpectralMatchIndexes;
         peptidesSearched = (results as TransientSearchEngineResults)!.PeptidesSearched;
         ReportProgress(new(100, "Finished Classic Search...", nestedIds));
     }

    private async Task<(TransientDatabaseContext Context, TransientDatabaseMetrics Metrics)> PerformPostSearchAnalysis(
        SpectralMatch[] allPsms,
        string outputFolder,
        List<string> nestedIds,
        string dbName,
        int totalProteins,
        HashSet<string> transientProteinAccessions,
        int transientPeptideCount,
        List<IBioPolymer> transientProteins,
        DbForTask transientDatabase,
        HashSet<int> updatedPsmIndexes)
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
        var transientPsms = FilterToTransientDatabaseOnly(psmList, transientProteinAccessions).ToList();
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
                CommonParameters,
                FileSpecificParameters,
                nestedIds);

            ProteinParsimonyResults proteinAnalysisResults =
                (ProteinParsimonyResults)await transientParsimonyEngine.RunAsync();
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
                        (ProteinScoringAndFdrResults)await new TransientProteinScoringAndFdrEngine(
                            baselineProteinGroups,
                            transientParsimonyGroups,
                            transientPsmsForParsimony,
                            _baselineProteinGroupFdr,
                            SearchParameters.NoOneHitWonders, SearchParameters.ModPeptidesAreDifferent,
                            true, CommonParameters, FileSpecificParameters, nestedIds).RunAsync();

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
                .Select(p => new CachedProteinGroup(p).CreateRuntimeCopy())
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
            CommonParameters = CommonParameters,
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

        return await Task.FromResult((analysisContext, aggregatedResult));
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
        await WritePsmsToTsvAsync(context.TransientPsms, transientPsmFile, SearchParameters.ModsToWriteSelection, false);
        FinishedWritingFile(transientPsmFile, nestedIds);

        if (writeAllResults)
        {
            string psmFile = Path.Combine(outputFolder,
                $"All{GlobalVariables.AnalyteType.GetSpectralMatchLabel()}s.{GlobalVariables.AnalyteType.GetSpectralMatchExtension()}");
            await WritePsmsToTsvAsync(context.AllPsms, psmFile, SearchParameters.ModsToWriteSelection, false);
            FinishedWritingFile(psmFile, nestedIds);
        }

        if (context.TransientPeptides is { Count: > 0 })
        {
            string transientPeptideFile = Path.Combine(outputFolder,
            $"{dbName}_All{GlobalVariables.AnalyteType}s.{GlobalVariables.AnalyteType.GetSpectralMatchExtension()}");
            await WritePsmsToTsvAsync(context.TransientPeptides, transientPeptideFile, SearchParameters.ModsToWriteSelection, true);
            FinishedWritingFile(transientPeptideFile, nestedIds);
        }

        if (writeAllResults)
        {
            string peptideFile = Path.Combine(outputFolder,
                $"All{GlobalVariables.AnalyteType}s.{GlobalVariables.AnalyteType.GetSpectralMatchExtension()}");
            await WritePsmsToTsvAsync(context.AllPeptides, peptideFile, SearchParameters.ModsToWriteSelection, true);
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

    private async Task WritePsmsToTsvAsync(IEnumerable<SpectralMatch> psms, string filePath, IReadOnlyDictionary<string, int> modstoWritePruned, bool writePeptideLevelResults = false)
    {
        await using StreamWriter output = new StreamWriter(filePath);
        bool includeOneOverK0Column = psms.Any(p => p.ScanOneOverK0.HasValue);
        await output.WriteLineAsync(SpectralMatch.GetTabSeparatedHeader(includeOneOverK0Column));
        foreach (var psm in psms)
        {
            await output.WriteLineAsync(psm.ToString(modstoWritePruned, writePeptideLevelResults, includeOneOverK0Column));
        }
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
        List<string> nestedIds, string taskId)
    {
        var transientDbList = new List<DbForTask> { transientDbPath };
        var transientDbLoader = new DatabaseLoadingEngine(CommonParameters,
            FileSpecificParameters, nestedIds, transientDbList, taskId,
            SearchParameters.DecoyType, SearchParameters.SearchTarget,
            LocalizableModificationTypes);
        var transientProteins = (transientDbLoader.Run() as DatabaseLoadingEngineResults)!.BioPolymers;
        return transientProteins;
    }

    /// <summary>
    /// Filters spectral matches to only include those that match exclusively to transient database proteins
    /// </summary>
    private IEnumerable<SpectralMatch> FilterToTransientDatabaseOnly(List<SpectralMatch> spectralMatches,
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
