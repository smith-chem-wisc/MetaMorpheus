#nullable enable
using EngineLayer.DiaLibrarySearch;

namespace TaskLayer;

/// <summary>
/// The settings a user sets for a DIA library search, as TOML. Only what a user should touch is here; the engine's
/// <see cref="DiaLibrarySearchParameters"/> (a positional record holding a run's calibration, which TOML cannot carry) is
/// built from these at run time, everything else at the engine's benchmarked defaults.
/// </summary>
public class DiaLibrarySearchTaskParameters
{
    /// <summary>Fragment m/z tolerance, ppm.</summary>
    public double FragmentTolerancePpm { get; set; } = new DiaLibrarySearchParameters().FragmentTolerancePpm;

    /// <summary>MS1 precursor m/z tolerance, ppm, around the run's calibrated MS1 offset.</summary>
    public double Ms1TolerancePpm { get; set; } = new DiaLibrarySearchParameters().Ms1TolerancePpm;

    /// <summary>Seed for the classifier, so a search can be repeated exactly.</summary>
    public int ClassifierSeed { get; set; } = new DiaLibrarySearchParameters().ClassifierSeed;

    /// <summary>The q-value at which precursors are counted in results.txt.</summary>
    public double QValueThreshold { get; set; } = 0.01;

    /// <summary>The engine's parameters for these settings, before a run's calibration is applied.</summary>
    public DiaLibrarySearchParameters ToEngineParameters() => new()
    {
        FragmentTolerancePpm = FragmentTolerancePpm,
        Ms1TolerancePpm = Ms1TolerancePpm,
        ClassifierSeed = ClassifierSeed,
    };
}
