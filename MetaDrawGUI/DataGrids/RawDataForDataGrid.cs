using System.IO;

namespace MetaDrawGUI
{
    class RawDataForDataGrid
    {
        #region Public Constructors

        public RawDataForDataGrid(string path)
        {
            FileName = Path.GetFileName(path);
            Use = true;
            FilePath = path;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Use { get; set; }
        public string FileName { get; private set; }
        //public string Parameters { get; private set; }
        //public bool InProgress { get; private set; }
        public string FilePath { get; private set; }

        #endregion Public Properties

        #region Public Methods

        /// <summary>
        /// Method to mark as in progress. Need the property setter to be private so user could not check off in GUI
        /// </summary>
        /// <param name="inProgress"></param>
        //public void SetInProgress(bool inProgress)
        //{
        //    InProgress = inProgress;
        //}

        //public void SetParametersText(string text)
        //{
        //    Parameters = text;
        //}

        #endregion Public Methods
    }
}
