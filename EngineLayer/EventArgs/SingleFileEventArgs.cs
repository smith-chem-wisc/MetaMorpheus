using System.Collections.Generic;

namespace EngineLayer
{
    public class SingleFileEventArgs : MyRecursiveEventArgs
    {
        public SingleFileEventArgs(string writtenFile, List<string> nestedIds) : base(nestedIds)
        {
            this.writtenFile = writtenFile;
        }

        public string writtenFile { get; private set; }
    }
}