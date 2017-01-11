using System;

namespace InternalLogic
{
    public class SingleFileEventArgs : EventArgs
    {
        public string writtenFile { get; private set; }

        public SingleFileEventArgs(string writtenFile)
        {
            this.writtenFile = writtenFile;
        }
    }
}