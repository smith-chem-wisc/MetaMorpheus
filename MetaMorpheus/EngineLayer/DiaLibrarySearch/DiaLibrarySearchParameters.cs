#nullable enable

namespace EngineLayer.DiaLibrarySearch;

/// <param name="FragmentTolerancePpm">How far an observed peak may sit from a library fragment and still be it.</param>
/// <param name="IrtHalfWindow">
/// How far, in iRT, a precursor's apex may fall from its library iRT. Wide enough here to absorb a run whose true
/// RT-to-iRT curve bows away from the linear map it was searched with.
/// </param>
public sealed record DiaLibrarySearchParameters(double FragmentTolerancePpm = 20, double IrtHalfWindow = 20);
