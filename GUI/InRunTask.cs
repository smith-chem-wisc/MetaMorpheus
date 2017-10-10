using TaskLayer;

namespace MetaMorpheusGUI
{
    public class InRunTask : ForTreeView
    {
        #region Public Fields

        public readonly MetaMorpheusTask task;

        #endregion Public Fields

        #region Public Constructors

        public InRunTask(string displayName, MetaMorpheusTask task) : base(displayName)
        {
            this.task = task;
        }

        #endregion Public Constructors
    }
}