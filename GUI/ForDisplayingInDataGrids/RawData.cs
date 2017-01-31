namespace MetaMorpheusGUI
{
    public class RawData
    {

        #region Public Constructors

        public RawData(string FileName)
        {
            this.FileName = FileName;
            if (FileName != null)
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