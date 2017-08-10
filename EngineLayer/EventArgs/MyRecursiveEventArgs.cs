using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class MyRecursiveEventArgs : EventArgs
    {
        #region Public Fields

        public readonly List<string> nestedIDs;

        #endregion Public Fields

        #region Public Constructors

        public MyRecursiveEventArgs(List<string> nestedIDs)
        {
            this.nestedIDs = nestedIDs;
        }

        #endregion Public Constructors
    }
}