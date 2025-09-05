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
using System.Windows;
using Easy.Common.Extensions;
using System.Windows.Data;

namespace GuiFunctions.MetaDraw;

public class MetaDrawDataLoader
{
    public event Action<string, int, int> ProgressChanged; // e.g. ("Loading Data Files", 2, 4)
    private readonly MetaDrawLogic _logic;
    private CancellationTokenSource? _cancellationTokenSource;

    public MetaDrawDataLoader(MetaDrawLogic logic)
    {
        _logic = logic;
        ProgressChanged += (stepName, current, total) =>
        {
            var step = LoadingProgressViewModel.Instance.Steps.FirstOrDefault(s => s.StepName == stepName);
            if (step != null)
            {
                step.Current = current;
                step.Total = total;

                // Only keep pop up open when base loading is still going on
                LoadingProgressViewModel.Instance.IsVisible = LoadingProgressViewModel.Instance.Steps.Any(s => s.Current < s.Total);
            }
        };
    }

    public async Task<List<string>> LoadAllAsync(
        bool loadSpectra,
        bool loadPsms,
        bool loadLibraries,
        ChimeraAnalysisTabViewModel? chimeraTabViewModel = null,
        BioPolymerTabViewModel? bioPolymerTabViewModel = null,
        DeconExplorationTabViewModel? deconExplorationTabViewModel = null)
    {
        // Cancel any previous loading operation
        _cancellationTokenSource?.Cancel();
        _cancellationTokenSource = new CancellationTokenSource();
        var cancellationToken = _cancellationTokenSource.Token;

        var allErrors = new List<string>();

        // Set up Progress Bars
        IProgress<(string, int, int)> progress = new Progress<(string, int, int)>(tuple => ProgressChanged?.Invoke(tuple.Item1, tuple.Item2, tuple.Item3));
        LoadingProgressViewModel.Instance.Steps.Clear();
        bool anySteps = false;
        if (loadSpectra && _logic.SpectraFilePaths.Count > 0)
        {
            LoadingProgressViewModel.Instance.Steps.Add(new LoadingStepViewModel { StepName = "Loading Spectra Files", Total = _logic.SpectraFilePaths.Count });
            anySteps = true;
        }
        if (loadPsms && _logic.SpectralMatchResultFilePaths.Count > 0)
        {
            LoadingProgressViewModel.Instance.Steps.Add(new LoadingStepViewModel { StepName = "Loading Search Results", Total = _logic.SpectralMatchResultFilePaths.Count });
            anySteps = true;
        }
        if (loadLibraries && _logic.SpectralLibraryPaths.Count > 0)
        {
            LoadingProgressViewModel.Instance.Steps.Add(new LoadingStepViewModel { StepName = "Loading Spectral Libraries", Total = _logic.SpectralLibraryPaths.Count });
            anySteps = true;
        }
        LoadingProgressViewModel.Instance.IsVisible = anySteps;

        if (chimeraTabViewModel is not null)
            chimeraTabViewModel.IsTabEnabled = false;
        if (bioPolymerTabViewModel is not null)
            bioPolymerTabViewModel.IsTabEnabled = false;
        if (deconExplorationTabViewModel is not null)
            deconExplorationTabViewModel.IsTabEnabled = false;

        // Load data in parallel
        var spectraTask = loadSpectra
            ? Task.Run(() => LoadSpectraFiles(progress, cancellationToken), cancellationToken)
            : Task.FromResult(new List<string>());
        var psmsTask = loadPsms
            ? Task.Run(() => LoadPsms(progress, haveLoadedSpectra: loadSpectra, cancellationToken), cancellationToken)
            : Task.FromResult(new List<string>());
        var librariesTask = loadLibraries
            ? Task.Run(() => LoadSpectralLibraries(progress, cancellationToken), cancellationToken)
            : Task.FromResult(new List<string>());

        List<string>[] results;
        try
        {
            results = await Task.WhenAll(spectraTask, psmsTask, librariesTask);
        }
        catch (OperationCanceledException)
        {
            allErrors.Add("Loading was canceled.");
            return allErrors;
        }

        if (cancellationToken.IsCancellationRequested)
        {
            allErrors.Add("Loading was canceled.");
            return allErrors;
        }

        _logic.FilterPsms();

        foreach (var errorList in results)
            allErrors.AddRange(errorList);

        // Load other tabs in the background. Only if there is data to load.
        if (bioPolymerTabViewModel != null && _logic.FilteredListOfPsms.Any())
        {
            _ = Task.Run(() =>
            {
                if (cancellationToken.IsCancellationRequested) 
                    return;

                if (!bioPolymerTabViewModel.IsDatabaseLoaded && !bioPolymerTabViewModel.DatabasePath.IsNullOrEmpty())
                    bioPolymerTabViewModel.LoadDatabaseCommand.Execute(null);

                if (bioPolymerTabViewModel.IsDatabaseLoaded)
                    bioPolymerTabViewModel.ProcessSpectralMatches(_logic.AllSpectralMatches, cancellationToken);

                bioPolymerTabViewModel.IsTabEnabled = true;
            }, cancellationToken);
        }
        if (deconExplorationTabViewModel != null && _logic.MsDataFiles.Any())
        {
            _ = Task.Run(() =>
            {
                if (cancellationToken.IsCancellationRequested) 
                    return;
                deconExplorationTabViewModel.MsDataFiles.Clear();
                for (int j = 0; j < _logic.MsDataFiles.Count; j++)
                {
                    deconExplorationTabViewModel.MsDataFiles.Add(_logic.MsDataFiles.ElementAt(j).Value);
                }

                deconExplorationTabViewModel.IsTabEnabled = true;
            }, cancellationToken);
        }
        if (chimeraTabViewModel != null && _logic.FilteredListOfPsms.Any() && _logic.MsDataFiles.Any())
        {
            _ = Task.Run(() =>
            {
                if (cancellationToken.IsCancellationRequested) 
                    return;

                chimeraTabViewModel.ProcessChimeraData(_logic.FilteredListOfPsms.ToList(), _logic.MsDataFiles, cancellationToken);

                chimeraTabViewModel.IsTabEnabled = true;
            }, cancellationToken);
        }

        return allErrors;
    }

    public List<string> LoadSpectraFiles(IProgress<(string, int, int)>? progress = null, CancellationToken cancellationToken = default)
    {
        var errors = new List<string>();
        int total = _logic.SpectraFilePaths.Count; 

        progress?.Report(("Loading Spectra Files", 0, total));
        for (int i = 0; i < total; i++)
        {
            if (cancellationToken.IsCancellationRequested)
                throw new OperationCanceledException();

            var filepath = _logic.SpectraFilePaths[i];
            var fileNameWithoutExtension = System.IO.Path.GetFileName(filepath.Replace(GlobalVariables.GetFileExtension(filepath), string.Empty));
            var spectraFile = MsDataFileReader.GetDataFile(filepath);
            if (spectraFile is TimsTofFileReader timsTofDataFile)
            {
                timsTofDataFile.LoadAllStaticData(maxThreads: Environment.ProcessorCount - 1);
            }
            else
            {
                spectraFile.InitiateDynamicConnection();
            }

            if (!_logic.MsDataFiles.TryAdd(fileNameWithoutExtension, spectraFile))
            {
                spectraFile.CloseDynamicConnection();
            }
            progress?.Report(("Loading Spectra Files", i + 1, total));
        }
        return errors;
    }

    private List<string> LoadPsms(IProgress<(string, int, int)>? progress = null, bool haveLoadedSpectra = false, CancellationToken cancellationToken = default)
    {
        var errors = new List<string>();
        int total = _logic.SpectralMatchResultFilePaths.Count;

        HashSet<string> fileNamesWithoutExtension = new HashSet<string>(
            _logic.SpectraFilePaths.Select(p => System.IO.Path.GetFileName(p.Replace(GlobalVariables.GetFileExtension(p), string.Empty))));
        List<SpectrumMatchFromTsv> psmsThatDontHaveMatchingSpectraFile = new List<SpectrumMatchFromTsv>();
        progress?.Report(("Loading Search Results", 0, total));
        for (int i = 0; i < total; i++)
        {
            if (cancellationToken.IsCancellationRequested)
                throw new OperationCanceledException();

            var resultsFile = _logic.SpectralMatchResultFilePaths[i];
            var psms = SpectrumMatchTsvReader.ReadTsv(resultsFile, out List<string> warnings);

            try
            {
                foreach (SpectrumMatchFromTsv psm in psms)
                {
                    if (fileNamesWithoutExtension.Contains(psm.FileNameWithoutExtension) || !haveLoadedSpectra)
                    {
                        _logic.AllSpectralMatches.Add(psm);
                    }
                    else
                    {
                        psmsThatDontHaveMatchingSpectraFile.Add(psm);
                    }

                    if (_logic.SpectralMatchesGroupedByFile.TryGetValue(psm.FileNameWithoutExtension, out var psmsForThisFile))
                    {
                        psmsForThisFile.Add(psm);
                    }
                    else
                    {
                        _logic.SpectralMatchesGroupedByFile.Add(psm.FileNameWithoutExtension, new ObservableCollection<SpectrumMatchFromTsv> { psm });
                    }
                }
            }
            catch (Exception e)
            {
                errors.Add("Error reading PSM file:\n" + e.Message);
            }
            progress?.Report(("Loading Search Results", i + 1, total));
        }

        if (psmsThatDontHaveMatchingSpectraFile.Any())
            foreach (var file in psmsThatDontHaveMatchingSpectraFile.GroupBy(p => p.FileNameWithoutExtension))
                errors.Add(file.Count() + " PSMs from " + file.Key + " were not loaded because this spectra file was not found");

        return errors;
    }

    private List<string> LoadSpectralLibraries(IProgress<(string, int, int)>? progress = null, CancellationToken cancellationToken = default)
    {
        var errors = new List<string>();
        int total = _logic.SpectralMatchResultFilePaths.Count; 

        progress?.Report(("Loading Spectral Libraries", 0, total));
        for (int i = 0; i < total; i++)
        {
            if (cancellationToken.IsCancellationRequested)
                throw new OperationCanceledException();

            try
            {
                _logic.SpectralLibrary = new SpectralLibrary(_logic.SpectralLibraryPaths.ToList());
            }
            catch (Exception e)
            {
                _logic.SpectralLibrary = null;
                errors.Add("Problem loading spectral library: " + e.Message);
            }
            progress?.Report(("Loading Spectral Libraries", i + 1, total));
        }
        return errors;
    }
}

