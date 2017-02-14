using System.Text;

namespace TaskLayer
{
    internal class MyGptmdTaskResults : MyTaskResults
    {

        #region Public Constructors

        public MyGptmdTaskResults(MetaMorpheusTask s) : base(s)
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            return sb.ToString();
        }

        #endregion Public Methods

    }
}