using System.Collections.Generic;

namespace EngineLayer
{
    public class SingleFileEventArgs : MyRecursiveEventArgs
    {

        #region Public Constructors

        public SingleFileEventArgs(string writtenFile, List<string> nestedIds) : base(nestedIds)
        {
            this.writtenFile = writtenFile;
        }

        #endregion Public Constructors

        #region Public Properties

        public string writtenFile { get; private set; }

        #endregion Public Properties

    }
}