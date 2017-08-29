using System.IO;

namespace MetaMorpheusGUI
{
    internal class RawDataForDataGrid
    {
        #region Public Constructors

        public RawDataForDataGrid(string fileName)
        {
            FileName = Path.GetFileNameWithoutExtension(fileName);
            Use = true;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Use { get; set; }
        public string FileName { get; }
        public bool InProgress { get; private set; }
        public string Parameters { get; set; }
        #endregion Public Properties

        #region Public Methods

        /// <summary>
        /// Method to mark as in progress. Need the property setter to be private so user could not check off in GUI
        /// </summary>
        /// <param name="inProgress"></param>
        public void SetInProgress(bool inProgress)
        {
            InProgress = inProgress;
        }

        #endregion Public Methods
    }
}