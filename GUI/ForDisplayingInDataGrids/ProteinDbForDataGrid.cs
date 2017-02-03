namespace MetaMorpheusGUI
{
    public class ProteinDbForDataGrid
    {

        #region Public Constructors

        public ProteinDbForDataGrid(string fileName)
        {
            Use = true;
            FileName = fileName;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Use { get; set; }
        public bool Contaminant { get; set; }
        public string FileName { get; private set; }
        public bool InProgress { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public void SetInProgress(bool inProgress)
        {
            InProgress = inProgress;
        }

        #endregion Public Methods

    }
}