using System;

namespace InternalLogicTaskLayer
{
    public class SingleTaskEventArgs : EventArgs
    {
        #region Public Constructors

        public SingleTaskEventArgs(MyTaskEngine theTask)
        {
            this.TheTask = theTask;
        }

        #endregion Public Constructors

        #region Public Properties

        public MyTaskEngine TheTask { get; private set; }

        #endregion Public Properties
    }
}