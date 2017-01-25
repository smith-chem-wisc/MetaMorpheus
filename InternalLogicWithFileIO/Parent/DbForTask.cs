namespace InternalLogicTaskLayer
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

        public string FileName { get; private set; }
        public bool IsContaminant { get; private set; }

        #endregion Public Properties

    }
}