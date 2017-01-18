using System;

namespace InternalLogicEngineLayer
{
    public class StringEventArgs : EventArgs
    {
        public string s { get; }
        public StringEventArgs(string s)
        {
            this.s = s;
        }
    }
}