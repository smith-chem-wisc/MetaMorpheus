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
            var sb = new StringBuilder();
            return sb.ToString();
        }
    }
}