using System.Text;

namespace InternalLogicTaskLayer
{
    internal class MyGPTMDTaskResults : MyTaskResults
    {
        public MyGPTMDTaskResults(MyTaskEngine s) : base(s)
        {
        }
        protected override string GetStringForOutput()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("MyGPTMDTaskResults:");
            return sb.ToString();
        }
    }
}