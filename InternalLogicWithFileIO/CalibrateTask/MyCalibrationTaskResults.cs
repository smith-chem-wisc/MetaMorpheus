using System.Text;

namespace InternalLogicTaskLayer
{
    internal class MyCalibrationTaskResults : MyTaskResults
    {
        public MyCalibrationTaskResults(MyTaskEngine s) : base(s)
        {
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("MyCalibrationTaskResults:");
            sb.AppendLine(base.ToString());
            return sb.ToString();
        }
    }
}