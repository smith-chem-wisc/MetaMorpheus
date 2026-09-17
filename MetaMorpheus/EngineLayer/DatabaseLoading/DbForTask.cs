#nullable enable
namespace EngineLayer.DatabaseLoading;

public class DbForTask
{
    public DbForTask(string filePath, bool isContaminant, string? decoyIdentifier = null)
    {
        FilePath = filePath;
        IsContaminant = isContaminant;
        FileName = System.IO.Path.GetFileName(filePath);
        var ext = GlobalVariables.GetFileExtension(filePath).ToLowerInvariant();
        IsSpectralLibrary = ext == ".msp" || ext == ".msl";
        DecoyIdentifier = decoyIdentifier ?? GlobalVariables.DecoyIdentifier;
    }

    /// <summary>
    /// The filename heuristic every front end uses to guess whether a database is a contaminant list.
    /// Ordinal, so it does not depend on the current culture. A guess only - callers should let the
    /// user override it.
    /// </summary>
    public static bool LooksLikeContaminant(string filePath) =>
        filePath.Contains("contaminant", System.StringComparison.OrdinalIgnoreCase)
        || filePath.Contains("crap", System.StringComparison.OrdinalIgnoreCase);

    public bool IsSpectralLibrary { get; }
    public string FilePath { get; }
    public bool IsContaminant { get; }
    public string FileName { get; }
    public string DecoyIdentifier { get; }
    public int? BioPolymerCount { get; internal set; } = null;
    public int? TargetCount { get; internal set; } = null;
    public int? DecoyCount { get; internal set; } = null;
}
