using EngineLayer;

namespace TaskLayer
{
    public interface IdbForTask
    {
        string FilePath { get; }
        bool IsContaminant { get; }
        string FileName { get; }
        bool IsSpectralLibrary { get; }
    }
    public class DbForTask
    {
        public DbForTask(string filePath, bool isContaminant)
        {
            FilePath = filePath;
            IsContaminant = isContaminant;
            FileName = System.IO.Path.GetFileName(filePath);
            IsSpectralLibrary = GlobalVariables.GetFileExtension(filePath).ToLowerInvariant() == ".msp";
        }
        public DbForTask(IdbForTask dbForTask)
        {
            this.dbForTask = dbForTask;
        }

        public bool IsSpectralLibrary { get; }
        public string FilePath { get; }
        public bool IsContaminant { get; }
        public string FileName { get; }

        private readonly IdbForTask dbForTask;
    }
}