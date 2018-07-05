using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class MyRecursiveEventArgs : EventArgs
    {
        public readonly List<string> nestedIDs;

        public MyRecursiveEventArgs(List<string> nestedIDs)
        {
            this.nestedIDs = nestedIDs;
        }
    }
}