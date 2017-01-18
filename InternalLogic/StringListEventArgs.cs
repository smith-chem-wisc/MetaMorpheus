using System;
using System.Collections.Generic;

namespace InternalLogicEngineLayer
{
    public class StringListEventArgs : EventArgs
    {
        public IEnumerable<string> StringList { get; }
        public StringListEventArgs(List<string> stringList)
        {
            this.StringList = stringList;
        }
    }
}