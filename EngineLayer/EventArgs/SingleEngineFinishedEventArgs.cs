using System;

namespace EngineLayer
{
    public class SingleEngineFinishedEventArgs : EventArgs
    {
        public readonly MetaMorpheusEngineResults myResults;

        public SingleEngineFinishedEventArgs(MetaMorpheusEngineResults myResults)
        {
            this.myResults = myResults;
        }

        public override string ToString()
        {
            return myResults.ToString();
        }
    }
}