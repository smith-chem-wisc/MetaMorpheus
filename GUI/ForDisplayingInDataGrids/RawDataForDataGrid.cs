namespace MetaMorpheusGUI
{
    internal class RawDataForDataGrid
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
        public string FileName { get; }
        public bool InProgress { get; set; }

        #endregion Public Properties

    }
}