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
            var sb = new StringBuilder();
            return sb.ToString();
        }
    }
}