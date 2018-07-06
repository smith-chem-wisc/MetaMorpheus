using System.IO;

namespace MetaMorpheusGUI
{
    public class RawDataForDataGrid
    {
        public RawDataForDataGrid(string path)
        {
            FileName = Path.GetFileName(path);
            Use = true;
            FilePath = path;
        }

        public bool Use { get; set; }
        public string FileName { get; private set; }
        public string Parameters { get; private set; }
        public bool InProgress { get; private set; }
        public string FilePath { get; private set; }

        /// <summary>
        /// Method to mark as in progress. Need the property setter to be private so user could not check off in GUI
        /// </summary>
        /// <param name="inProgress"></param>
        public void SetInProgress(bool inProgress)
        {
            InProgress = inProgress;
        }

        public void SetParametersText(string text)
        {
            Parameters = text;
        }
    }
}