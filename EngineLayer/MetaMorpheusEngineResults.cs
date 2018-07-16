using System;
using System.Text;

namespace EngineLayer
{
    public class MetaMorpheusEngineResults
    {
        internal TimeSpan Time;

        public MetaMorpheusEngineResults(MetaMorpheusEngine s)
        {
            MyEngine = s;
        }

        public MetaMorpheusEngine MyEngine { get; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append("Time to run: " + Time);
            return sb.ToString();
        }
    }
}