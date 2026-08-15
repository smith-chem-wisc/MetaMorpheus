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
internal static class ToolkitNamespaceIsNotShadowed
{
    internal static Avalonia.Point Origin => new(0, 0);
    internal static Avalonia.Media.Color Black => Avalonia.Media.Colors.Black;
}
