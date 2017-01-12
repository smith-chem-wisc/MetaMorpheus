using System;

namespace InternalLogicTaskLayer
{
    public class SingleTaskEventArgs : EventArgs
    {
        #region Public Constructors

        public SingleTaskEventArgs(MyTaskEngine theTask)
        {
            this.theTask = theTask;
        }

        #endregion Public Constructors

        #region Public Properties

        public MyTaskEngine theTask { get; private set; }

        #endregion Public Properties
    }
}