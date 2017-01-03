using System;
using System.Text;

namespace IndexSearchAndAnalyze
{
    public abstract class MyResults
    {
        protected MyParams s;

        public MyResults(MyParams s)
        {
            this.s = s;
        }

        public TimeSpan Time { get; internal set; }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("Results:");
            sb.AppendLine();
            sb.Append("Time to run: " + Time);
            return sb.ToString();
        }
    }
}