using System.Collections.Generic;

namespace EngineLayer
{
    public class StringEventArgs : MyRecursiveEventArgs
    {
        public StringEventArgs(string s, List<string> nestedIDs)
            : base(nestedIDs)
        {
            this.S = s;
        }

        public string S { get; }
    }
}