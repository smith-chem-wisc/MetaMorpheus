using System;
using System.Collections.Generic;

namespace InternalLogicEngineLayer
{
    public class StringListEventArgs : EventArgs
    {

        #region Public Constructors

        public StringListEventArgs(List<string> stringList)
        {
            this.StringList = stringList;
        }

        #endregion Public Constructors

        #region Public Properties

        public IEnumerable<string> StringList { get; }

        #endregion Public Properties

    }
}