using System.Text;

namespace InternalLogicTaskLayer
{
    internal class MyCalibrationTaskResults : MyTaskResults
    {

        #region Public Constructors

        public MyCalibrationTaskResults(MyTaskEngine s) : base(s)
        {
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override string GetStringForOutput()
        {
            var sb = new StringBuilder();
            return sb.ToString();
        }

        #endregion Protected Methods

    }
}