using System;

namespace InternalLogic
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