using System;

namespace EngineLayer
{
    public class SingleEngineEventArgs : EventArgs
    {
        public SingleEngineEventArgs(MetaMorpheusEngine myEngine)
        {
            MyEngine = myEngine;
        }

        public MetaMorpheusEngine MyEngine { get; private set; }
    }
}