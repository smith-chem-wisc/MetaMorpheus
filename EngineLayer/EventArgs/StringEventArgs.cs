using System;

namespace EngineLayer
{
    public class StringEventArgs : EventArgs
    {

        #region Public Constructors

        public StringEventArgs(string s)
        {
            this.s = s;
        }

        #endregion Public Constructors

        #region Public Properties

        public string s { get; }

        #endregion Public Properties

    }
}