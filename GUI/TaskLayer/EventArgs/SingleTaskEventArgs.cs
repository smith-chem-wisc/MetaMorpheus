using System;

namespace TaskLayer
{
    public class SingleTaskEventArgs : EventArgs
    {
        #region Public Constructors

        public SingleTaskEventArgs(string displayName)
        {
            this.DisplayName = displayName;
        }

        #endregion Public Constructors

        #region Public Properties

        public string DisplayName { get; private set; }

        #endregion Public Properties
    }
}