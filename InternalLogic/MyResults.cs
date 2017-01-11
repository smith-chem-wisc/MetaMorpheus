using System;
using System.Text;

namespace InternalLogicEngineLayer
{
    public abstract class MyResults
    {
        protected MyEngine s { get; private set; }

        public MyResults(MyEngine s)
        {
            this.s = s;
        }

        public TimeSpan Time { get; internal set; }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("Time to run: " + Time);
            return sb.ToString();
        }
    }
}