using TaskLayer;

namespace MetaMorpheusGUI
{
    internal class ProteinDbForDataGrid
    {

        #region Public Constructors

        public ProteinDbForDataGrid(string fileName)
        {
            Use = true;
            FileName = fileName;
            if (fileName.ToUpper().Contains("contaminant".ToUpper())
                || fileName.ToUpper().Contains("crap".ToUpper()))
                Contaminant = true;
        }

        public ProteinDbForDataGrid(DbForTask uu)
        {
            Use = true;
            Contaminant = uu.IsContaminant;
            FileName = uu.FileName;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool Use { get; set; }
        public bool Contaminant { get; set; }
        public string FileName { get; }
        public bool InProgress { get; set; }

        #endregion Public Properties

    }
}