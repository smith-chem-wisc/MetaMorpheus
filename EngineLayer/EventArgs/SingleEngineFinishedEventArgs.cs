using System;

namespace EngineLayer
{
    public class SingleEngineFinishedEventArgs : EventArgs
    {
        public readonly MetaMorpheusEngineResults MyResults;

        public SingleEngineFinishedEventArgs(MetaMorpheusEngineResults myResults)
        {
            MyResults = myResults;
        }

        public override string ToString()
        {
            return MyResults.ToString();
        }
    }
}