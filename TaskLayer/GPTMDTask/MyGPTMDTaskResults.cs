using System.Text;

namespace TaskLayer
{
    internal class MyGPTMDTaskResults : MyTaskResults
    {

        #region Public Constructors

        public MyGPTMDTaskResults(MetaMorpheusTask s) : base(s)
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