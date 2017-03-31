using System;

namespace EngineLayer
{
    public class SingleEngineEventArgs : EventArgs
    {

        #region Public Constructors

        public SingleEngineEventArgs(MetaMorpheusEngine myEngine)
        {
            this.myEngine = myEngine;
        }

        #endregion Public Constructors

        #region Public Properties

        public MetaMorpheusEngine myEngine { get; private set; }

        #endregion Public Properties

    }
}