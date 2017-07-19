namespace TaskLayer
{
    public class DbForTask
    {
        #region Public Constructors

        public DbForTask(string fileName, bool isContaminant)
        {
            FilePath = fileName;
            IsContaminant = isContaminant;
        }

        #endregion Public Constructors

        #region Public Properties

        public string FilePath { get; }
        public bool IsContaminant { get; }

        #endregion Public Properties
    }
}