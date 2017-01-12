using System;

namespace InternalLogicEngineLayer
{
    public class SingleFileEventArgs : EventArgs
    {
        #region Public Constructors

        public SingleFileEventArgs(string writtenFile)
        {
            this.writtenFile = writtenFile;
        }

        #endregion Public Constructors

        #region Public Properties

        public string writtenFile { get; private set; }

        #endregion Public Properties
    }
}