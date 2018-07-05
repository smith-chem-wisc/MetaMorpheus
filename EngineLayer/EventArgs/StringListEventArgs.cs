using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class StringListEventArgs : EventArgs
    {
        public StringListEventArgs(List<string> stringList)
        {
            this.StringList = stringList;
        }

        public IEnumerable<string> StringList { get; }
    }
}