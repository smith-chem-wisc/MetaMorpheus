using System;

namespace InternalLogicWithFileIO
{
    public class SingleTaskEventArgs : EventArgs
    {
        public MyTaskEngine theTask { get; private set; }

        public SingleTaskEventArgs(MyTaskEngine theTask)
        {
            this.theTask = theTask;
        }
    }
}