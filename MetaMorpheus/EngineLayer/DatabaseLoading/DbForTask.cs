#nullable enable
namespace EngineLayer.DatabaseLoading;

public class DbForTask
{
    public DbForTask(string filePath, bool isContaminant, string? decoyIdentifier = null)
    {
        FilePath = filePath;
        IsContaminant = isContaminant;
        FileName = System.IO.Path.GetFileName(filePath);
        IsSpectralLibrary = GlobalVariables.GetFileExtension(filePath).ToLowerInvariant() == ".msp";
        DecoyIdentifier = decoyIdentifier ?? GlobalVariables.DecoyIdentifier;
    }

    public bool IsSpectralLibrary { get; }
    public string FilePath { get; }
    public bool IsContaminant { get; }
    public string FileName { get; }
    public string DecoyIdentifier { get; }
    public int? BioPolymerCount { get; internal set; } = null;
    public int? TargetCount { get; internal set; } = null;
    public int? DecoyCount { get; internal set; } = null;
}
