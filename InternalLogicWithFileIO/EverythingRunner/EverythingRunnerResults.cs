using System;
using InternalLogicEngineLayer;
using System.Text;

namespace InternalLogicTaskLayer
{
    internal class EverythingRunnerResults : MyResults
    {
        public EverythingRunnerResults(MyEngine s) : base(s)
        {
        }

        protected override string GetStringForOutput()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("EverythingRunnerResults:");
            return sb.ToString();
        }
    }
}