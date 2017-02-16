using TaskLayer;

namespace MetaMorpheusGUI
{
    internal class MetaMorpheusTaskForDataGrid
    {

        #region Public Fields

        public readonly MetaMorpheusTask metaMorpheusTask;

        #endregion Public Fields

        #region Public Constructors

        public MetaMorpheusTaskForDataGrid(MetaMorpheusTask theTask)
        {
            this.metaMorpheusTask = theTask;
        }

        #endregion Public Constructors

        #region Public Properties

        public string Task { get { return metaMorpheusTask.GetType().Name; } }
        public bool InProgress { get; set; }

        #endregion Public Properties

    }
}