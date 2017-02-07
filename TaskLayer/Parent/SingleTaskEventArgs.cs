using System;

namespace TaskLayer
{
    public class SingleTaskEventArgs : EventArgs
    {
        #region Public Constructors

        public SingleTaskEventArgs(MetaMorpheusTask theTask)
        {
            this.TheTask = theTask;
        }

        #endregion Public Constructors

        #region Public Properties

        public MetaMorpheusTask TheTask { get; private set; }

        #endregion Public Properties
    }
}