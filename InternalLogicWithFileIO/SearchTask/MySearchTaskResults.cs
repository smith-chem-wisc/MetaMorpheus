using InternalLogicEngineLayer;
using System.Text;

namespace InternalLogicTaskLayer
{
    internal class MySearchTaskResults : MyTaskResults
    {
        public MySearchTaskResults(MyTaskEngine s) : base(s)
        {
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("MySearchTaskResults:");
            sb.Append(base.ToString());
            return sb.ToString();
        }
    }
}