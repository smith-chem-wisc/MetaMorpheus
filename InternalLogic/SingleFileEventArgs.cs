using System;

namespace InternalLogicEngineLayer
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