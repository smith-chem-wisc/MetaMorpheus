namespace TaskLayer
{
    public class DbForTask
    {
        #region Public Constructors

        public DbForTask(string fileName, bool isContaminant)
        {
            FileName = fileName;
            IsContaminant = isContaminant;
        }

        #endregion Public Constructors

        #region Public Properties

        public string FileName { get; }
        public bool IsContaminant { get; }

        #endregion Public Properties
    }
}