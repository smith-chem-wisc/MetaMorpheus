using System;

namespace InternalLogicEngineLayer
{
    public class SingleEngineEventArgs : EventArgs
    {
        private MyEngine myEngine;

        public SingleEngineEventArgs(MyEngine myEngine)
        {
            this.myEngine = myEngine;
        }
    }
}