using System;

namespace MetaMorpheusLogic
{
    public class SingleTaskEventArgs : EventArgs
    {
        public MyTask theTask { get; private set; }

        public SingleTaskEventArgs(MyTask theTask)
        {
            this.theTask = theTask;
        }
    }
}