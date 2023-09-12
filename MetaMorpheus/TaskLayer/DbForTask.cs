using EngineLayer;

namespace TaskLayer
{
    public class DbForTask
    {
        public DbForTask(string filePath, bool isContaminant)
        {
            FilePath = filePath;
            IsContaminant = isContaminant;
            FileName = System.IO.Path.GetFileName(filePath);
            IsSpectralLibrary = GlobalVariables.GetFileExtension(filePath).ToLowerInvariant() == ".msp";
        }

        public bool IsSpectralLibrary { get; }
        public string FilePath { get; }
        public bool IsContaminant { get; }
        public string FileName { get; }
    }
}