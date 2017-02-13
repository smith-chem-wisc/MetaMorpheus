using System.Text;

namespace TaskLayer
{
    internal class MyCalibrationTaskResults : MyTaskResults
    {

        #region Public Constructors

        public MyCalibrationTaskResults(MetaMorpheusTask s) : base(s)
        {
        }

        #endregion Public Constructors

        #region Protected Properties
        
        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            return sb.ToString();
        }

        #endregion Protected Properties

    }
}