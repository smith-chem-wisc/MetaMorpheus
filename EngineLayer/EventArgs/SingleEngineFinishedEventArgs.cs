using System;

namespace EngineLayer
{
    public class SingleEngineFinishedEventArgs : EventArgs
    {
        #region Public Fields

        public readonly MetaMorpheusEngineResults myResults;

        #endregion Public Fields

        #region Public Constructors

        public SingleEngineFinishedEventArgs(MetaMorpheusEngineResults myResults)
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