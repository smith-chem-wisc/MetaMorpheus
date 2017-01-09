using System;

namespace IndexSearchAndAnalyze
{
    public class ProgressEventArgs : EventArgs

    {
        public int new_progress;
        public string v;

        public ProgressEventArgs(int new_progress, string v)
        {
            this.new_progress = new_progress;
            this.v = v;
        }
    }
}