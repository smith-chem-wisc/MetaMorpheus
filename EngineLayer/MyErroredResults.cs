using System.Text;

namespace EngineLayer
{
    public class MyErroredResults : MetaMorpheusEngineResults
    {

        #region Private Fields

        private readonly string errorString;

        #endregion Private Fields

        #region Public Constructors

        public MyErroredResults(MetaMorpheusEngine s, string errorString) : base(s)
        {
            this.errorString = errorString;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append(errorString);
            return sb.ToString();
        }

        #endregion Public Methods

    }
}