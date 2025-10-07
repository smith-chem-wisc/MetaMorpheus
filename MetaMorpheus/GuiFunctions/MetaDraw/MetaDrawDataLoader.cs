﻿#nullable enable
using EngineLayer;
using Readers.SpectralLibrary;
using Readers;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Data;
using System.Collections.Concurrent;
using System.IO;
using Readers.InternalResults;
using MassSpectrometry;
using Nett;
using TaskLayer;

namespace GuiFunctions.MetaDraw;

public class MetaDrawDataLoader
{
    private readonly MetaDrawLogic _logic;
    private CancellationTokenSource? _cancellationTokenSource;

    public MetaDrawDataLoader(MetaDrawLogic logic)
    {
        _logic = logic;
        BindingOperations.EnableCollectionSynchronization(LoadingProgressViewModel.Instance.Steps, _progressLock);
    }

    public async Task<List<string>> LoadAllAsync(
        bool loadSpectra,
        bool loadPsms,
        bool loadLibraries,
        ChimeraAnalysisTabViewModel? chimeraTabViewModel = null,
        BioPolymerTabViewModel? bioPolymerTabViewModel = null,
        DeconExplorationTabViewModel? deconExplorationTabViewModel = null)
    {
        // Cancel any previous run
        _cancellationTokenSource?.Cancel();
        _cancellationTokenSource = new CancellationTokenSource();
        var token = _cancellationTokenSource.Token;
        var allErrors = new List<string>();

        // Prepare progress visualization
        PrepareProgressVisualization(loadSpectra, loadPsms, loadLibraries, _logic);
        if (chimeraTabViewModel is not null) 
            chimeraTabViewModel.IsTabEnabled = false;
        if (bioPolymerTabViewModel is not null) 
            bioPolymerTabViewModel.IsTabEnabled = false;
        if (deconExplorationTabViewModel is not null) 
            deconExplorationTabViewModel.IsTabEnabled = false;

        // Start core artifact loads (independent steps run concurrently)
        var spectraTask = loadSpectra 
            ? LoadSpectraAsync(token) 
            : Task.FromResult(new List<string>());
        var psmsTask = loadPsms 
            ? LoadPsmsAsync(loadSpectra, token) 
            : Task.FromResult(new List<string>());
        var librariesTask = loadLibraries 
            ? LoadLibrariesAsync(token) 
            : Task.FromResult(new List<string>());
        var proseAndTomlTask = _logic.SpectralMatchResultFilePaths.Any()
            ? TryLoadProseAndSearchToml(bioPolymerTabViewModel, deconExplorationTabViewModel)
            : Task.FromResult((false, false));

        List<string>[] results;
        try 
        { 
            results = await Task.WhenAll(spectraTask, psmsTask, librariesTask).ConfigureAwait(false); 
        }
        catch (OperationCanceledException) { return new() { "Loading was canceled." }; }

        allErrors.AddRange(results.SelectMany(e => e));
        if (token.IsCancellationRequested) return allErrors;

        string outputDirectory = null!;
        if (_logic.SpectralMatchResultFilePaths.Count > 0)
        {
            outputDirectory = Path.GetDirectoryName(_logic.SpectralMatchResultFilePaths.First()) ?? string.Empty;
        }
        else if (_logic.SpectraFilePaths.Count > 0)
        {
            outputDirectory = Path.GetDirectoryName(_logic.SpectraFilePaths.First()) ?? string.Empty;
        }
        else
            outputDirectory = Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);

        await proseAndTomlTask;
        // Tabs: run in background, post minimal changes to UI
        if (bioPolymerTabViewModel is not null && _logic.AllSpectralMatches.Any())
        {
            Task.Run(() => ProcessBioPolymerTab(bioPolymerTabViewModel, outputDirectory, token), token);
        }
        if (deconExplorationTabViewModel is not null && _logic.MsDataFiles.Any())
        {
            Task.Run(() => ProcessDeconTab(deconExplorationTabViewModel, outputDirectory, token), token);
        }
        if (chimeraTabViewModel is not null && _logic.FilteredListOfPsms.Any() && _logic.MsDataFiles.Any())
        {
            Task.Run(() => ProcessChimeraTab(chimeraTabViewModel, outputDirectory, token), token);
        }

        return allErrors;
    }

    #region Core loaders (non-blocking)

    public async Task<List<string>> LoadSpectraAsync(CancellationToken token)
    {
        var errors = new ConcurrentBag<string>();
        var files = _logic.SpectraFilePaths.Where(path => File.Exists(path) || Directory.Exists(path)).ToArray();
        var total = files.Length;
        ReportProgress(StepSpectra, 0, total);

        var maxDeg = Math.Max(2, Math.Min(Environment.ProcessorCount, 8));
        var completed = 0;

        await Task.Run(() =>
        {
            Parallel.ForEach(files,
                new ParallelOptions { MaxDegreeOfParallelism = maxDeg, CancellationToken = token },
                path =>
                {
                    token.ThrowIfCancellationRequested();

                    var ext = GlobalVariables.GetFileExtension(path);
                    var key = System.IO.Path.GetFileName(path.Replace(ext, string.Empty));

                    var dataFile = MsDataFileReader.GetDataFile(path);
                    try
                    {
                        if (dataFile is TimsTofFileReader tims)
                        {
                            tims.LoadAllStaticData(maxThreads: Math.Max(1, Environment.ProcessorCount - 1));
                        }
                        else
                        {
                            dataFile.InitiateDynamicConnection();
                        }

                        lock (_logic.ThreadLocker) // MsDataFiles is a Dictionary
                        {
                            if (!_logic.MsDataFiles.ContainsKey(key))
                            {
                                _logic.MsDataFiles.Add(key, dataFile);
                            }
                            else
                            {
                                dataFile.CloseDynamicConnection();
                            }
                        }
                    }
                    catch (Exception ex)
                    {
                        try { dataFile.CloseDynamicConnection(); } catch { /* best effort */ }
                        errors.Add($"Error loading spectra file '{path}': {ex.Message}");
                    }
                    finally
                    {
                        var done = Interlocked.Increment(ref completed);
                        ReportProgress(StepSpectra, done, total);
                    }
                });
        }, token).ConfigureAwait(false);

        return errors.ToList();
    }

    public async Task<List<string>> LoadPsmsAsync(bool requireMatchingSpectra, CancellationToken token)
    {
        var errors = new ConcurrentBag<string>();
        var tsvs = _logic.SpectralMatchResultFilePaths;
        var total = tsvs.Count;
        ReportProgress(StepPsms, 0, total);

        var spectraNames = new HashSet<string>(
            _logic.SpectraFilePaths.Select(p =>
            {
                var ext = GlobalVariables.GetFileExtension(p);
                return System.IO.Path.GetFileName(p.Replace(ext, string.Empty));
            }),
            StringComparer.OrdinalIgnoreCase);

        var acceptedAll = new ConcurrentBag<SpectrumMatchFromTsv>();
        var acceptedByFile = new ConcurrentDictionary<string, List<SpectrumMatchFromTsv>>();
        var missingDueToNoSpectra = new ConcurrentBag<SpectrumMatchFromTsv>();
        var maxDeg = Math.Max(2, Math.Min(Environment.ProcessorCount, 8));
        var completed = 0;

        await Task.Run(() =>
        {
            Parallel.ForEach(tsvs,
                new ParallelOptions { MaxDegreeOfParallelism = maxDeg, CancellationToken = token },
                rf =>
                {
                    token.ThrowIfCancellationRequested();
                    try
                    {
                        var psms = SpectrumMatchTsvReader.ReadTsv(rf, out var warnings);
                        foreach (var w in warnings)
                        {
                            errors.Add($"Warning parsing '{rf}': {w}");
                        }

                        foreach (var p in psms)
                        {
                            bool hasMatchingSpectra = spectraNames.Contains(p.FileNameWithoutExtension);
                            if (!requireMatchingSpectra || hasMatchingSpectra)
                            {
                                acceptedAll.Add(p);
                                acceptedByFile.AddOrUpdate(p.FileNameWithoutExtension,
                                    _ => new List<SpectrumMatchFromTsv> { p },
                                    (_, list) => { list.Add(p); return list; });
                            }
                            else
                            {
                                missingDueToNoSpectra.Add(p);
                            }
                        }
                    }
                    catch (Exception ex) { errors.Add($"Error reading PSM file:\n{ex}"); }
                    finally
                    {
                        var done = Interlocked.Increment(ref completed);
                        ReportProgress(StepPsms, done, total);
                    }
                });
        }, token).ConfigureAwait(false);

        // Merge once; FilteredListOfPsms & grouped dict are already collection-synchronized
        lock (_logic.ThreadLocker)
        {
            _logic.AllSpectralMatches.AddRange(acceptedAll
                .OrderByDescending(p => p is PsmFromTsv { BetaPeptideScore: not null } psm 
                    ? p.Score + psm.BetaPeptideScore.Value
                    : p.Score
                )
                .ThenByDescending(p => p.DeltaScore));
            foreach (var kv in acceptedByFile)
            {
                var ordered = kv.Value
                    .OrderByDescending(p => p is PsmFromTsv { BetaPeptideScore: not null } psm
                        ? p.Score + psm.BetaPeptideScore.Value
                        : p.Score
                    )
                    .ThenByDescending(p => p.DeltaScore);
                if (_logic.SpectralMatchesGroupedByFile.TryGetValue(kv.Key, out var existing))
                {
                    foreach (var p in ordered) 
                        existing.Add(p);
                }
                else
                {
                    _logic.SpectralMatchesGroupedByFile[kv.Key] = new(ordered);
                }
            }

            _logic.FilterPsms(); 
        }

        // Missing spectra messages (after merge)
        var missingMsgs = missingDueToNoSpectra
            .GroupBy(p => p.FileNameWithoutExtension, StringComparer.OrdinalIgnoreCase)
            .Select(g => $"{g.Count()} PSMs from {g.Key} were not loaded because this spectra file was not found");

        foreach (var m in missingMsgs) 
            errors.Add(m);
        return errors.ToList();
    }

    public async Task<List<string>> LoadLibrariesAsync(CancellationToken token)
    {
        var errors = new List<string>();
        var total = Math.Max(1, _logic.SpectralLibraryPaths.Count);
        ReportProgress(StepLibraries, 0, total);

        try
        {
            await Task.Run(() =>
            {
                token.ThrowIfCancellationRequested();
                if (_logic.SpectralLibraryPaths.Count > 0)
                {
                    _logic.SpectralLibrary = new SpectralLibrary(_logic.SpectralLibraryPaths.ToList());
                }
            }, token).ConfigureAwait(false);
        }
        catch (OperationCanceledException) { throw; }
        catch (Exception e)
        {
            _logic.SpectralLibrary = null;
            errors.Add("Problem loading spectral library: " + e.Message);
        }
        finally
        {
            ReportProgress(StepLibraries, total, total);
        }

        return errors;
    }

    public async Task<(bool Database, bool SearchParams)> TryLoadProseAndSearchToml(BioPolymerTabViewModel? bpTabVm, DeconExplorationTabViewModel? deconTabVm)
    {
        bool loadedDb = false;
        bool loadedSearchParams = false;

        var searchDirectories = _logic.SpectralMatchResultFilePaths
            .Select(Path.GetDirectoryName)
            .Where(p => p != null)
            .Select(p =>
            {
                // If directory name contains "Individual File Results", go up one parent
                if (Path.GetFileName(p)?.Contains("Individual File Results", StringComparison.OrdinalIgnoreCase) ?? false)
                {
                    var parent = Directory.GetParent(p!);
                    return parent?.FullName ?? p;
                }
                return p;
            })
            .Distinct()
            .ToArray();

        // Try to find search tomls to load decon parameters. 
        try
        {
            HashSet<string?> distinctSearchTomls = new();
            foreach (var searchDir in searchDirectories)
            {
                if (searchDir == null) 
                    continue;

                var parentDir = Directory.GetParent(searchDir)?.FullName ?? null;
                var searchDirName = Path.GetFileName(searchDir);

                if (parentDir == null) 
                    continue;

                // Attempt to locate search toml
                var settingsDir = Directory.GetDirectories(parentDir).FirstOrDefault(p => p.Contains("Task Settings"));
                if (settingsDir != null)
                {
                    distinctSearchTomls.Add(
                        Directory.GetFiles(settingsDir, $"*{searchDirName}config.toml").FirstOrDefault());
                }
            }

            List<DeconvolutionParameters> precursorParameters = new();
            List<DeconvolutionParameters> productParameters = new();

            foreach (var searchToml in distinctSearchTomls.Where(p => p != null))
            {
                var task = Toml.ReadFile<SearchTask>(searchToml, MetaMorpheusTask.tomlConfig);
                precursorParameters.Add(task.CommonParameters.PrecursorDeconvolutionParameters);
                productParameters.Add(task.CommonParameters.ProductDeconvolutionParameters);
            }

            // IF we find params, and have multiple, take the last. 
            // Not the best way to do it, but without an equality comparison on decon params, it's hard to know if they are different or identical. 
            if (precursorParameters.Any() && productParameters.Any())
            {
                var precursor = precursorParameters.Last();
                var product = productParameters.Last();
                MetaDrawSettingsViewModel.Instance.DeconHostViewModel = new(precursor, product);
                loadedSearchParams = true;

                if (deconTabVm != null) 
                {
                    deconTabVm.DeconHostViewModel = new(precursor, product);
                }
            }
        }
        catch (Exception) { loadedSearchParams = false; }

        if (bpTabVm is null)
            return (loadedDb, loadedSearchParams);

        string[] acceptableDbTypes = new[] { ".fasta", ".fa", ".xml" };

        // Try to find databases from the prose files. 
        try
        {
            List<string> databases = new();
            foreach (var searchDir in searchDirectories)
            {
                if (searchDir == null) 
                    continue;
                var proseFile = MetaMorpheusProseFile.LocateInDirectory(searchDir);
                if (proseFile == null) 
                    continue;
                databases.AddRange(proseFile.DatabasePaths);
            }

            // Db tab only supports one database at a time currently, pick the one which occurs the most, then first in order
            foreach (var dbGroup in databases
                .Where(p => File.Exists(p) && acceptableDbTypes.Any(p.Contains))
                            .GroupBy(db => db, StringComparer.OrdinalIgnoreCase)
                            .OrderByDescending(g => g.Count())
                            .ThenBy(g => databases.IndexOf(g.Key)))
            {
                bpTabVm.DatabasePath = dbGroup.Key;
                loadedDb = true;
                break;
            }
        }
        catch (Exception) { loadedDb = false; }


        return (loadedDb, loadedSearchParams);
    }

    #endregion

    #region Tab Processing 

    private void ProcessBioPolymerTab(BioPolymerTabViewModel tab, string outDir, CancellationToken token)
    {
        if (token.IsCancellationRequested) return;

        // Load DB if present but not loaded
        if (!tab.IsDatabaseLoaded && !string.IsNullOrEmpty(tab.DatabasePath))
        {
            tab.LoadDatabaseCommand.Execute(null);
        }

        if (tab.IsDatabaseLoaded)
        {
            tab.ProcessSpectralMatches(_logic.AllSpectralMatches, token);
        }

        tab.ExportDirectory = outDir;
        tab.IsTabEnabled = true;
    }

    private void ProcessDeconTab(DeconExplorationTabViewModel tab, string outDir, CancellationToken token)
    {
        if (token.IsCancellationRequested) 
            return;

        tab.MsDataFiles.Clear();
        foreach (var file in _logic.MsDataFiles.Values)
        {
            tab.MsDataFiles.Add(file);
        }

        tab.ExportDirectory = outDir;
        tab.IsTabEnabled = true;
    }

    private void ProcessChimeraTab(ChimeraAnalysisTabViewModel tab, string outDir, CancellationToken token)
    {
        if (token.IsCancellationRequested) 
            return;

        tab.ProcessChimeraData(_logic.FilteredListOfPsms.ToList(), _logic.MsDataFiles, token);

        tab.ExportDirectory = outDir;
        tab.IsTabEnabled = true;
    }

    #endregion

    #region UI Progress

    private readonly object _progressLock = new();
    private readonly Dictionary<string, LoadingStepViewModel> _steps = new(StringComparer.Ordinal);
    private const string StepSpectra = "Loading Spectra Files";
    private const string StepPsms = "Loading Search Results";
    private const string StepLibraries = "Loading Spectral Libraries";

    private void PrepareProgressVisualization(
        bool loadSpectra, bool loadPsms, bool loadLibraries, MetaDrawLogic logic)
    {
        var vm = LoadingProgressViewModel.Instance;
        lock (_progressLock)
        {
            vm.Steps.Clear();
            _steps.Clear();

            if (loadSpectra && logic.SpectraFilePaths.Count > 0)
                AddStep(StepSpectra, logic.SpectraFilePaths.Count);

            if (loadPsms && logic.SpectralMatchResultFilePaths.Count > 0)
                AddStep(StepPsms, logic.SpectralMatchResultFilePaths.Count);

            if (loadLibraries && logic.SpectralLibraryPaths.Count > 0)
                AddStep(StepLibraries, Math.Max(1, logic.SpectralLibraryPaths.Count));

            vm.IsVisible = vm.Steps.Any();
        }
    }

    private void AddStep(string name, int total)
    {
        var step = new LoadingStepViewModel { StepName = name, Total = total, Current = 0 };
        LoadingProgressViewModel.Instance.Steps.Add(step);
        _steps[name] = step;
    }

    private void ReportProgress(string stepName, int current, int total)
    {
        lock (_progressLock)
        {
            if (_steps.TryGetValue(stepName, out var s))
            {
                s.Total = total;
                s.Current = current;
            }
            LoadingProgressViewModel.Instance.IsVisible =
                LoadingProgressViewModel.Instance.Steps.Any(st => st.Current < st.Total);
        }
    }

    #endregion
}

