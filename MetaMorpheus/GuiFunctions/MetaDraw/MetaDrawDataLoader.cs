#nullable enable
using EngineLayer;
using MassSpectrometry;
using Readers.SpectralLibrary;
using Readers;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TaskLayer;

namespace GuiFunctions.MetaDraw;

public class MetaDrawDataLoader
{
    public event Action<string, int, int> ProgressChanged; // e.g. ("Loading Data Files", 2, 4)
    public event Action<string, bool> TabReadyChanged; // e.g. ("Chimera", true)
    private readonly MetaDrawLogic _logic;

    public MetaDrawDataLoader(MetaDrawLogic logic)
    {
        _logic = logic;
    }
    public async Task<List<string>> LoadAllAsync(
        bool loadSpectra,
        bool loadPsms,
        bool loadLibraries,
        ChimeraAnalysisTabViewModel? chimeraTabViewModel = null,
        BioPolymerTabViewModel? bioPolymerTabViewModel = null,
        DeconExplorationTabViewModel? deconExplorationTabViewModel = null)
    {
        var allErrors = new List<string>();
        var progress = new Progress<(string, int, int)>(tuple => ProgressChanged?.Invoke(tuple.Item1, tuple.Item2, tuple.Item3));

        var tasks = new List<Task<List<string>>>();

        // Load Base Information
        if (loadSpectra)
            tasks.Add(LoadSpectraFilesAsync(progress));
        if (loadPsms)
            tasks.Add(LoadPsmsAsync(progress, loadSpectra));
        if (loadLibraries)
            tasks.Add(LoadSpectralLibrariesAsync(progress));
        
        // Aggregate errors
        var results = await Task.WhenAll(tasks);
        foreach (var errorList in results)
            allErrors.AddRange(errorList);

        // Post-load processing
        if (chimeraTabViewModel != null)
        {
            chimeraTabViewModel.ProcessChimeraData(_logic.FilteredListOfPsms.ToList(), _logic.MsDataFiles);
            TabReadyChanged?.Invoke("Chimera", true);
        }
        if (bioPolymerTabViewModel != null && bioPolymerTabViewModel.IsDatabaseLoaded) // TODO: Load database if available
        {
            bioPolymerTabViewModel.ProcessSpectralMatches(_logic.AllSpectralMatches);
            TabReadyChanged?.Invoke("BioPolymer", true);
        }
        if (deconExplorationTabViewModel != null)
        {
            deconExplorationTabViewModel.MsDataFiles.Clear();
            foreach (var dataFile in _logic.MsDataFiles)
            {
                deconExplorationTabViewModel.MsDataFiles.Add(dataFile.Value);
            }
            TabReadyChanged?.Invoke("Deconvolution", true);
        }

        return allErrors;
    }

    public async Task<List<string>> LoadSpectraFilesAsync(IProgress<(string, int, int)> progress = null)
    {
        var errors = new List<string>();
        int total = _logic.SpectraFilePaths.Count;
        progress?.Report(("Loading Data Files", 0, total));
        for (int i = 0; i < total; i++)
        {
            var filepath = _logic.SpectraFilePaths[i];
            var fileNameWithoutExtension = System.IO.Path.GetFileName(filepath.Replace(GlobalVariables.GetFileExtension(filepath), string.Empty));
            var spectraFile = MsDataFileReader.GetDataFile(filepath);
            if (spectraFile is TimsTofFileReader timsTofDataFile)
            {
                timsTofDataFile.LoadAllStaticData(maxThreads: Environment.ProcessorCount - 1); // timsTof files need to load all static data before they can be used, as dynamic access is not available for them
            }
            else
            {
                spectraFile.InitiateDynamicConnection();
            }

            if (!_logic.MsDataFiles.TryAdd(fileNameWithoutExtension, spectraFile))
            {
                spectraFile.CloseDynamicConnection();
                // print warning? but probably unnecessary. this means the data file was loaded twice. 
                // which is an error but not an important one because the data is loaded
            }
            progress?.Report(("Loading Data Files", i + 1, total));
        }
        return errors;
    }

    public async Task<List<string>> LoadPsmsAsync(IProgress<(string, int, int)> progress = null, bool haveLoadedSpectra = false)
    {
        var errors = new List<string>();
        int total = _logic.SpectralMatchResultFilePaths.Count;

        HashSet<string> fileNamesWithoutExtension = new HashSet<string>(
            _logic.SpectraFilePaths.Select(p => System.IO.Path.GetFileName(p.Replace(GlobalVariables.GetFileExtension(p), string.Empty))));
        List<SpectrumMatchFromTsv> psmsThatDontHaveMatchingSpectraFile = new List<SpectrumMatchFromTsv>();
        progress?.Report(("Loading Search Results", 0, total));
        for (int i = 0; i < total; i++)
        {
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

        _logic.FilterPsms();
        return errors;
    }

    public async Task<List<string>> LoadSpectralLibrariesAsync(IProgress<(string, int, int)> progress = null)
    {
        var errors = new List<string>();
        int total = _logic.SpectralMatchResultFilePaths.Count;
        progress?.Report(("Loading Search Results", 0, total));
        for (int i = 0; i < total; i++)
        {
            try
            {
                _logic.SpectralLibrary = new SpectralLibrary(_logic.SpectralLibraryPaths.ToList());
            }
            catch (Exception e)
            {
                _logic.SpectralLibrary = null;
                errors.Add("Problem loading spectral library: " + e.Message);
            }
            progress?.Report(("Loading Search Results", i + 1, total));
        }
        return errors;
    }
}

