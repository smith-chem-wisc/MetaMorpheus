using System.IO;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using GuiFunctions;
using TaskLayer;

namespace MetaMorpheusGUI
{
    public class ProteinDbForDataGrid : BaseViewModel
    {
        #region Public Constructors

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
            DecoyIdentifier = uu.DecoyIdentifier;
        }

        #endregion Public Constructors

        #region Public Properties

        private bool _use;
        private bool _isContaminant;
        private bool _inProgress;
        private string _decoyIdentifier = GlobalVariables.DecoyIdentifier;

        public bool Use
        {
            get => _use;
            set { _use = value; OnPropertyChanged(nameof(Use)); }
        }
        public bool Contaminant
        {
            get => _isContaminant;
            set { _isContaminant = value; OnPropertyChanged(nameof(Contaminant)); }
        }
        public string FileName { get; private set; }
        public string FilePath { get; private set; }
        public bool InProgress
        {
            get => _inProgress;
            private set { _inProgress = value; OnPropertyChanged(nameof(InProgress)); }
        }

        public string DecoyIdentifier
        {
            get => _decoyIdentifier;
            set { _decoyIdentifier = value; OnPropertyChanged(nameof(DecoyIdentifier)); }
        }

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