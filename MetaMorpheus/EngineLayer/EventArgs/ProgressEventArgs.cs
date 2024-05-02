using System;
using System.Collections.Generic;

namespace EngineLayer
{
    [Serializable]
    public class ProgressEventArgs : MyRecursiveEventArgs

    {
        public int NewProgress;
        public string V;

        public ProgressEventArgs(int newProgress, string v, List<string> nestedIDs) : base(nestedIDs)
        {
            NewProgress = newProgress;
            V = v;
        }
    }
}