using System.Text;

namespace InternalLogicTaskLayer
{
    internal class MySearchTaskResults : MyTaskResults
    {
        public MySearchTaskResults(MyTaskEngine s) : base(s)
        {
        }

        protected override string StringForOutput
        {
            get
            {
                var sb = new StringBuilder();
                return sb.ToString();
            }
        }
    }
}