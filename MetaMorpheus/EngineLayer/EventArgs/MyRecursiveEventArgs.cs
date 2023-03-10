using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class MyRecursiveEventArgs : EventArgs
    {
        public readonly List<string> NestedIDs;

        public MyRecursiveEventArgs(List<string> nestedIDs)
        {
            this.NestedIDs = nestedIDs;
        }
    }
}