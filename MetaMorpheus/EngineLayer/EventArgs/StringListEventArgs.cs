using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class StringListEventArgs : EventArgs
    {
        public StringListEventArgs(List<string> stringList)
        {
            StringList = stringList;
        }

        public IEnumerable<string> StringList { get; }
    }
}