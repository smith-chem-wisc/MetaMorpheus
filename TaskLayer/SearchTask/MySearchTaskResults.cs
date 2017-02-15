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

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(base.ToString());
            return sb.ToString();
        }

        #endregion Public Methods

    }
}