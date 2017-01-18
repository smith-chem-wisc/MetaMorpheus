using System;

namespace InternalLogicEngineLayer
{
    public class ProgressEventArgs : EventArgs

    {
        public int new_progress;
        public string v;

        public ProgressEventArgs(int newProgress, string v)
        {
            this.new_progress = newProgress;
            this.v = v;
        }
    }
}