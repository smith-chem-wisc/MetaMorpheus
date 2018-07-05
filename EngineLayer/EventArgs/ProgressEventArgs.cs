using System.Collections.Generic;

namespace EngineLayer
{
    public class ProgressEventArgs : MyRecursiveEventArgs

    {
        public int new_progress;
        public string v;

        public ProgressEventArgs(int newProgress, string v, List<string> nestedIDs) : base(nestedIDs)
        {
            this.new_progress = newProgress;
            this.v = v;
        }
    }
}