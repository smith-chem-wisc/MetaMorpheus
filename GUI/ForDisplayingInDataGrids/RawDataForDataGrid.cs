namespace MetaMorpheusGUI
{
    class RawDataForDataGrid
    {

        #region Public Constructors

        public RawDataForDataGrid(string fileName)
        {
            FileName = fileName;
            Use = true;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Use { get; set; }
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