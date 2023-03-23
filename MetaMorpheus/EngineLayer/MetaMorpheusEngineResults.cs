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
            sb.Append("Engine type: " + MyEngine.GetType().Name + "\n");
            sb.Append("Time to run engine: " + Time);
            return sb.ToString();
        }
    }
}