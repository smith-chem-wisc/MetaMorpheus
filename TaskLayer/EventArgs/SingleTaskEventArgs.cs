using System;

namespace TaskLayer
{
    public class SingleTaskEventArgs : EventArgs
    {

        #region Public Constructors

        public SingleTaskEventArgs(string theTask)
        {
            this.TaskId = theTask;
        }

        #endregion Public Constructors

        #region Public Properties

        public string TaskId { get; private set; }

        #endregion Public Properties

    }
}