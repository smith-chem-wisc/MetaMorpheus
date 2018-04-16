namespace TaskLayer
{
    public class DbForTask
    {
        #region Public Constructors

        public DbForTask(string filePath, bool isContaminant)
        {
            FilePath = filePath;
            IsContaminant = isContaminant;
            FileName = System.IO.Path.GetFileName(filePath);
        }

        #endregion Public Constructors

        #region Public Properties

        public string FilePath { get; }
        public bool IsContaminant { get; }
        public string FileName { get; }

        #endregion Public Properties
    }
}