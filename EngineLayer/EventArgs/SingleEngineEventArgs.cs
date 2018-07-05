using System;

namespace EngineLayer
{
    public class SingleEngineEventArgs : EventArgs
    {
        public SingleEngineEventArgs(MetaMorpheusEngine myEngine)
        {
            this.myEngine = myEngine;
        }

        public MetaMorpheusEngine myEngine { get; private set; }
    }
}