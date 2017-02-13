using System.Text;

namespace TaskLayer
{
    internal class MySearchTaskResults : MyTaskResults
    {

        #region Public Constructors

        public MySearchTaskResults(MetaMorpheusTask s) : base(s)
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