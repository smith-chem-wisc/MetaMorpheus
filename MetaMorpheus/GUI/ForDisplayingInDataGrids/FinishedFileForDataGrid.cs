namespace MetaMorpheusGUI
{
    public class FinishedFileForDataGrid
    {
        #region Public Constructors

        public FinishedFileForDataGrid(string filePath)
        {
            FilePath = filePath;
        }

        #endregion Public Constructors

        #region Public Properties

        public string FilePath { get; private set; }

        #endregion Public Properties
    }
}