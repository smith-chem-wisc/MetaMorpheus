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

        protected override string StringForOutput
        {
            get
            {
                var sb = new StringBuilder();
                return sb.ToString();
            }
        }

        #endregion Protected Properties

    }
}