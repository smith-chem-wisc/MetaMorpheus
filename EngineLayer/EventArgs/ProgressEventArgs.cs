using System;

namespace EngineLayer
{
    public class ProgressEventArgs : EventArgs

    {

        #region Public Fields

        public int new_progress;
        public string v;

        #endregion Public Fields

        #region Public Constructors

        public ProgressEventArgs(int newProgress, string v)
        {
            this.new_progress = newProgress;
            this.v = v;
        }

        #endregion Public Constructors

    }
}