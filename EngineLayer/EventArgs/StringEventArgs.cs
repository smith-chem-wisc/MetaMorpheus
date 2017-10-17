using System.Collections.Generic;

namespace EngineLayer
{
    public class StringEventArgs : MyRecursiveEventArgs
    {
        #region Public Constructors

        public StringEventArgs(string s, List<string> nestedIDs)
            : base(nestedIDs)
        {
            this.S = s;
        }

        #endregion Public Constructors

        #region Public Properties

        public string S { get; }

        #endregion Public Properties
    }
}