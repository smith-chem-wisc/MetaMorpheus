using System;
using System.Collections.Generic;
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
    public async Task<List<string>> LoadAllAsync(bool loadSpectra, bool loadPsms, bool loadLibraries, bool loadChimera, bool loadBioPolymer)
    {
        var allErrors = new List<string>();
        var progress = new Progress<(string, int, int)>(tuple => ProgressChanged?.Invoke(tuple.Item1, tuple.Item2, tuple.Item3));

        var tasks = new List<Task<List<string>>>();

        if (loadSpectra)
            tasks.Add(_logic.LoadSpectraFilesAsync(progress));
        if (loadPsms)
            tasks.Add(_logic.LoadPsmsAsync(progress));
        if (loadLibraries)
            tasks.Add(_logic.LoadSpectralLibrariesAsync(progress));

        // Chimera and BioPolymer can be loaded in parallel if needed
        // You may want to trigger TabReadyChanged when each finishes

        var results = await Task.WhenAll(tasks);
        foreach (var errorList in results)
            allErrors.AddRange(errorList);

        // Signal tab readiness
        if (loadSpectra) TabReadyChanged?.Invoke("Deconvolution", true);
        if (loadPsms) TabReadyChanged?.Invoke("BioPolymer", true);
        if (loadChimera) TabReadyChanged?.Invoke("Chimera", true);

        return allErrors;
    }
}

