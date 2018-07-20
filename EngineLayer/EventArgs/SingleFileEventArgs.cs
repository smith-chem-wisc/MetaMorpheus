using System.Collections.Generic;

namespace EngineLayer
{
    public class SingleFileEventArgs : MyRecursiveEventArgs
    {
        public SingleFileEventArgs(string writtenFile, List<string> nestedIds) : base(nestedIds)
        {
            WrittenFile = writtenFile;
        }

        public string WrittenFile { get; private set; }
    }
}