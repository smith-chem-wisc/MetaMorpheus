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
            if (s.Level <= 1)
            {
                sb.AppendLine(GetType().Name.ToString() + ":");
                sb.AppendLine(GetStringForOutput());
                sb.AppendLine("\tTime to run: " + Time);
            }
            else
            {
                sb.AppendLine("\t" + GetType().Name.ToString() + ":");
                sb.AppendLine(GetStringForOutput());
                sb.Append("\t\tTime to run: " + Time);
            }
            return sb.ToString();
        }

        protected abstract string GetStringForOutput();
    }
}