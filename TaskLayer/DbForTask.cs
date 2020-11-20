namespace TaskLayer
{
    public class DbForTask
    {
        public DbForTask(string filePath, bool isContaminant)
        {
            FilePath = filePath;
            IsContaminant = isContaminant;
            FileName = System.IO.Path.GetFileName(filePath);
        }

        public DbForTask(string filePath, bool isContaminant, bool isSpectralLibrary)
        {
            FilePath = filePath;
            IsContaminant = isContaminant;
            IsSpectralLibrary = isSpectralLibrary;
            FileName = System.IO.Path.GetFileName(filePath);
        }

        public bool IsSpectralLibrary { get; }
        public string FilePath { get; }
        public bool IsContaminant { get; }
        public string FileName { get; }
    }
}