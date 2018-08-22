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

        public string FilePath { get; }
        public bool IsContaminant { get; }
        public string FileName { get; }
    }
}