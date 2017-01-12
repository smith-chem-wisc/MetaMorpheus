using System.Text;

namespace InternalLogicTaskLayer
{
    internal class MyGPTMDTaskResults : MyTaskResults
    {
        public MyGPTMDTaskResults(MyTaskEngine s) : base(s)
        {
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("MyGPTMDTaskResults:");
            sb.Append(base.ToString());
            return sb.ToString();
        }
    }
}