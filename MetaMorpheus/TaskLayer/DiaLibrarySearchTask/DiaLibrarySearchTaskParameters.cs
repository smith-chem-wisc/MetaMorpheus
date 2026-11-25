#nullable enable
using System.Collections.Generic;
using EngineLayer.DiaLibrarySearch;

namespace TaskLayer;

/// <summary>How a predicted library's decoys are made (<see cref="DiaLibrarySearchTaskParameters.LibraryDecoys"/>).</summary>
public enum DiaLibraryDecoys
{
    /// <summary>A composition-preserving permutation of each target protein's pieces, named DECOY_&lt;accession&gt;_f1.</summary>
    Permuted,

    /// <summary>Each target protein reversed, named DECOY_&lt;accession&gt;.</summary>
    Reversed,
}

/// <summary>
/// The settings a user sets for a DIA library search, as TOML. Only what a user should touch is here; the engine's
/// <see cref="DiaLibrarySearchParameters"/> (a positional record holding a run's calibration, which TOML cannot carry) is
/// built from these at run time, everything else at the engine's benchmarked defaults.
/// </summary>
public class DiaLibrarySearchTaskParameters
{
    /// <summary>
    /// Fragment m/z tolerance, ppm. 20 by default, for Orbitrap MS2; set 10 for Orbitrap Astral MS2, which gained 2.5-13% on
    /// PXD071205's 5-min runs (see <see cref="DiaLibrarySearchParameters.FragmentTolerancePpm"/>).
    /// </summary>
    public double FragmentTolerancePpm { get; set; } = new DiaLibrarySearchParameters().FragmentTolerancePpm;

    /// <summary>MS1 precursor m/z tolerance, ppm, around the run's calibrated MS1 offset.</summary>
    public double Ms1TolerancePpm { get; set; } = new DiaLibrarySearchParameters().Ms1TolerancePpm;

    /// <summary>Seed for the classifier, so a search can be repeated exactly.</summary>
    public int ClassifierSeed { get; set; } = new DiaLibrarySearchParameters().ClassifierSeed;

    /// <summary>
    /// The precursor charges a predicted library holds, when the search is given proteins and no library
    /// (<see cref="DiaPredictedLibrary"/>). 2 and 3 by default, as in the benchmark libraries; Prosit takes 1 to 6.
    /// </summary>
    public List<int> LibraryPrecursorCharges { get; set; } = [2, 3];

    /// <summary>
    /// The normalized collision energy a predicted library's spectra are predicted at: set the runs' HCD energy. 27 by
    /// default, the energy of the HeLa benchmark runs (PXD005573, PXD022589); PXD071205's Astral runs were 25.
    /// </summary>
    public int LibraryNce { get; set; } = 27;

    /// <summary>
    /// How a predicted library's decoys are made. Permuted by default: each decoy protein a composition-preserving
    /// permutation of its target's tryptic pieces (mzLib's entrapment generator), so every target peptide has a near-isomeric
    /// decoy, as in the benchmark libraries. Reversed decoys are for a digestion the permuter refuses (not fully specific,
    /// cleavage-blocking modifications respected, or a rule such as trypsin|P).
    /// </summary>
    public DiaLibraryDecoys LibraryDecoys { get; set; } = DiaLibraryDecoys.Permuted;

    /// <summary>The engine's parameters for these settings, before a run's calibration is applied.</summary>
    public DiaLibrarySearchParameters ToEngineParameters() => new()
    {
        FragmentTolerancePpm = FragmentTolerancePpm,
        Ms1TolerancePpm = Ms1TolerancePpm,
        ClassifierSeed = ClassifierSeed,
    };
}
