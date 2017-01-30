using System;

namespace EngineLayer
{
    public class SingleEngineEventArgs : EventArgs
    {
        #region Public Constructors

        public SingleEngineEventArgs(MyEngine myEngine)
        {
            this.myEngine = myEngine;
        }

        #endregion Public Constructors

        #region Public Properties

        public MyEngine myEngine { get; private set; }

        #endregion Public Properties
    }
}