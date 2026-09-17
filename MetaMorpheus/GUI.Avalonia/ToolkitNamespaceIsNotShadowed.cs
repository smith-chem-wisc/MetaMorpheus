namespace MetaMorpheusAvalonia;

/// <summary>
/// A compile-time guard, not runtime code. These are fully-qualified toolkit references with no
/// using directive, which only resolve while this assembly's root namespace does not contain a
/// member called Avalonia.
///
/// The root namespace used to be MetaMorpheus.Avalonia, and inside it "Avalonia.Point" bound to
/// MetaMorpheus.Avalonia.Point and failed with CS0234 naming the wrong namespace, fixable only with
/// global::. Nothing hit it because every file had a using directive above the namespace - but the
/// MetaDraw port brings ~600 references to Avalonia types, and it would have hit it constantly.
/// </summary>
// Excluded from coverage deliberately, matching how the WPF GUI is treated.
//
// PRECEDENT: Test.csproj references CMD, EngineLayer, GuiFunctions and TaskLayer - but NOT
// GUI.csproj. The WPF views are therefore never loaded by the test host, never instrumented, and
// contribute nothing to the coverage denominator: Codecov reports 0 files under MetaMorpheus/GUI/
// and 79 under MetaMorpheus/GuiFunctions/. That split is the whole point of GuiFunctions - the
// logic was moved out of the views so that it COULD be covered, leaving behind only the parts that
// cannot meaningfully be unit tested.
//
// GUI.Avalonia cannot use the same mechanism, because its tests construct real windows headlessly
// and so must reference this project. [ExcludeFromCodeCoverage] is the equivalent, and is already
// the house instrument for untestable code (~20 files on master carry it, e.g. the design-time
// view models in GuiFunctions/ViewModels/Deconvolution).
//
// THE RULE THIS ENCODES: if something here can be tested, it does not belong here - it belongs on
// the view model, where it is covered like the rest of GuiFunctions. That is why SpectraPatterns
// and DatabasePatterns now live on MainWindowViewModel rather than in this file. Do not add logic
// to a class carrying this attribute in order to avoid writing a test for it.
// This one is a compile-time guard with no runtime callers by design, so it could never be covered.
[System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
internal static class ToolkitNamespaceIsNotShadowed
{
    internal static Avalonia.Point Origin => new(0, 0);
    internal static Avalonia.Media.Color Black => Avalonia.Media.Colors.Black;
}
