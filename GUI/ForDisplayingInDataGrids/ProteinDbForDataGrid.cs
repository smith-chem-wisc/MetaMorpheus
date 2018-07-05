using System.IO;
using TaskLayer;

namespace MetaMorpheusGUI
{
    internal class ProteinDbForDataGrid
    {
        public ProteinDbForDataGrid(string FilePath)
        {
            Use = true;
            this.FilePath = FilePath;
            if (FilePath.ToUpperInvariant().Contains("contaminant".ToUpperInvariant())
                || FilePath.ToUpperInvariant().Contains("CRAP"))
            {
                Contaminant = true;
            }
            FileName = Path.GetFileName(FilePath);
        }

        public ProteinDbForDataGrid(DbForTask uu)
        {
            Use = true;
            Contaminant = uu.IsContaminant;
            FilePath = uu.FilePath;
            FileName = uu.FileName;
        }

        public bool Use { get; set; }
        public bool Contaminant { get; set; }
        public string FileName { get; private set; }
        public string FilePath { get; private set; }
        public bool InProgress { get; private set; }

        /// <summary>
        /// Method to mark as in progress. Need the property setter to be private so user could not check off in GUI
        /// </summary>
        /// <param name="inProgress"></param>
        public void SetInProgress(bool inProgress)
        {
            InProgress = inProgress;
        }
    }
}