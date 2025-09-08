#nullable enable
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
        if (chimeraTabViewModel is not null) chimeraTabViewModel.IsTabEnabled = false;
        if (bioPolymerTabViewModel is not null) bioPolymerTabViewModel.IsTabEnabled = false;
        if (deconExplorationTabViewModel is not null) deconExplorationTabViewModel.IsTabEnabled = false;

        // Start core artifact loads (independent steps run concurrently)
        var spectraTask = loadSpectra ? LoadSpectraAsync(token) : Task.FromResult(new List<string>());
        var psmsTask = loadPsms ? LoadPsmsAsync(loadSpectra, token) : Task.FromResult(new List<string>());
        var librariesTask = loadLibraries ? LoadLibrariesAsync(token) : Task.FromResult(new List<string>());

        List<string>[] results;
        try { results = await Task.WhenAll(spectraTask, psmsTask, librariesTask).ConfigureAwait(false); }
        catch (OperationCanceledException) { return new() { "Loading was canceled." }; }

        allErrors.AddRange(results.SelectMany(e => e));
        if (token.IsCancellationRequested) return allErrors;

        string outputDirectory = null!;
        if (_logic.SpectralMatchResultFilePaths.Count > 0)
            outputDirectory = Path.GetDirectoryName(_logic.SpectralMatchResultFilePaths.First()) ?? string.Empty;
        else if (_logic.SpectraFilePaths.Count > 0)
            outputDirectory = Path.GetDirectoryName(_logic.SpectraFilePaths.First()) ?? string.Empty;
        else
            outputDirectory = Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments);

        // Tabs: run in background, post minimal changes to UI
        if (bioPolymerTabViewModel is not null && _logic.AllSpectralMatches.Any())
            Task.Run(() => ProcessBioPolymerTab(bioPolymerTabViewModel, outputDirectory, token), token);
        if (deconExplorationTabViewModel is not null && _logic.MsDataFiles.Any())
            Task.Run(() => ProcessDeconTab(deconExplorationTabViewModel, outputDirectory, token), token);
        if (chimeraTabViewModel is not null && _logic.FilteredListOfPsms.Any() && _logic.MsDataFiles.Any())
            Task.Run(() => ProcessChimeraTab(chimeraTabViewModel, outputDirectory, token), token);

        return allErrors;
    }


    #region Core loaders (non-blocking)

    private async Task<List<string>> LoadSpectraAsync(CancellationToken token)
    {
        var errors = new ConcurrentBag<string>();
        var files = _logic.SpectraFilePaths;
        var total = files.Count;
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
                            tims.LoadAllStaticData(maxThreads: Math.Max(1, Environment.ProcessorCount - 1));
                        else
                            dataFile.InitiateDynamicConnection();

                        lock (_logic.ThreadLocker) // MsDataFiles is a Dictionary
                        {
                            if (!_logic.MsDataFiles.ContainsKey(key))
                                _logic.MsDataFiles.Add(key, dataFile);
                            else
                                dataFile.CloseDynamicConnection();
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

    private async Task<List<string>> LoadPsmsAsync(bool requireMatchingSpectra, CancellationToken token)
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
                        foreach (var w in warnings) errors.Add($"Warning parsing '{rf}': {w}");

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
                    catch (Exception ex)
                    {
                        errors.Add("Error reading PSM file:\n" + ex.Message);
                    }
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
            _logic.AllSpectralMatches.AddRange(acceptedAll);
            foreach (var kv in acceptedByFile)
            {
                if (_logic.SpectralMatchesGroupedByFile.TryGetValue(kv.Key, out var existing))
                {
                    foreach (var p in kv.Value) 
                        existing.Add(p);
                }
                else
                {
                    _logic.SpectralMatchesGroupedByFile[kv.Key] = new ObservableCollection<SpectrumMatchFromTsv>(kv.Value);
                }
            }

            _logic.FilterPsms(); 
        }

        // Missing spectra messages (after merge)
        var missingMsgs = missingDueToNoSpectra
            .GroupBy(p => p.FileNameWithoutExtension, StringComparer.OrdinalIgnoreCase)
            .Select(g => $"{g.Count()} PSMs from {g.Key} were not loaded because this spectra file was not found");

        foreach (var m in missingMsgs) errors.Add(m);
        return errors.ToList();
    }

    private async Task<List<string>> LoadLibrariesAsync(CancellationToken token)
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
                    _logic.SpectralLibrary = new SpectralLibrary(_logic.SpectralLibraryPaths.ToList());
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

    #endregion

    #region Tab Processing 

    private void ProcessBioPolymerTab(BioPolymerTabViewModel tab, string outDir, CancellationToken token)
    {
        if (token.IsCancellationRequested) return;

        // Load DB if present but not loaded
        if (!tab.IsDatabaseLoaded && !string.IsNullOrEmpty(tab.DatabasePath))
            tab.LoadDatabaseCommand.Execute(null);

        if (tab.IsDatabaseLoaded)
            tab.ProcessSpectralMatches(_logic.AllSpectralMatches, token);

        tab.ExportDirectory = outDir;
        tab.IsTabEnabled = true;
    }

    private void ProcessDeconTab(DeconExplorationTabViewModel tab, string outDir, CancellationToken token)
    {
        if (token.IsCancellationRequested) 
            return;

        tab.MsDataFiles.Clear();
        foreach (var file in _logic.MsDataFiles.Values)
            tab.MsDataFiles.Add(file);

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

