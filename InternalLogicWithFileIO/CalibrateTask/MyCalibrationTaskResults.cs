using System.Text;

namespace InternalLogicTaskLayer
{
    internal class MyCalibrationTaskResults : MyTaskResults
    {
        public MyCalibrationTaskResults(MyTaskEngine s) : base(s)
        {
        }
        protected override string GetStringForOutput()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("MyCalibrationTaskResults:");
            return sb.ToString();
        }
    }
}