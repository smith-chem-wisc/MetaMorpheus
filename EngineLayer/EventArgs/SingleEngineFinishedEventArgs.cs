using System;

namespace EngineLayer
{
    public class SingleEngineFinishedEventArgs : EventArgs
    {

        #region Private Fields

        private readonly MyResults myResults;

        #endregion Private Fields

        #region Public Constructors

        public SingleEngineFinishedEventArgs(MyResults myResults)
        {
            this.myResults = myResults;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return myResults.ToString();
        }

        #endregion Public Methods

    }
}