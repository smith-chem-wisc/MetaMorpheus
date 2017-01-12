using InternalLogicEngineLayer;
using System.Text;

namespace InternalLogicTaskLayer
{
    internal class MySearchTaskResults : MyTaskResults
    {
        public MySearchTaskResults(MyTaskEngine s) : base(s)
        {
        }
        protected override string GetStringForOutput()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("MySearchTaskResults:");
            return sb.ToString();
        }
    }
}