using TaskLayer;

namespace MetaMorpheusGUI
{
    internal class PreRunTask
    {
        #region Public Fields

        public readonly MetaMorpheusTask metaMorpheusTask;

        #endregion Public Fields

        #region Public Constructors

        public PreRunTask(MetaMorpheusTask theTask)
        {
            this.metaMorpheusTask = theTask;
        }

        #endregion Public Constructors

        #region Public Properties

        public string DisplayName { get; set; }

        #endregion Public Properties

        #region Public Methods

        public PreRunTask Clone()
        {
            return (PreRunTask)this.MemberwiseClone();
        }

        #endregion Public Methods
    }
}