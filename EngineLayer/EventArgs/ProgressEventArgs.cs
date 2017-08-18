using System.Collections.Generic;

namespace EngineLayer
{
    public class ProgressEventArgs : MyRecursiveEventArgs

    {
        #region Public Fields

        public int new_progress;
        public string v;

        #endregion Public Fields

        #region Public Constructors

        public ProgressEventArgs(int newProgress, string v, List<string> nestedIDs) : base(nestedIDs)
        {
            this.new_progress = newProgress;
            this.v = v;
        }

        #endregion Public Constructors
    }
}