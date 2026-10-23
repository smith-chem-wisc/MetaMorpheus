#nullable enable

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// One library precursor scored against one DIA run: its apex, where that falls in iRT, and its score. Deliberately not a
/// <see cref="SpectralMatch"/>: a DIA identification is a chromatographic peak group, not one scan's hypothesis.
/// </summary>
/// <param name="PrecursorIndex">Index of the precursor in the library searched.</param>
/// <param name="QValue">
/// Target-decoy q-value among targets; NaN for a decoy, which competes but is not itself reported.
/// </param>
public sealed record DiaPrecursorMatch(
    int PrecursorIndex,
    string FullSequence,
    int Charge,
    double PrecursorMz,
    bool IsDecoy,
    Irt LibraryIrt,
    RtMinutes ApexRt,
    Irt ApexIrt,
    double Score,
    double QValue = double.NaN);
