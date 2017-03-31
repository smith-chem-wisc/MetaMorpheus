using System.Text;

namespace EngineLayer
{
    public class MyErroredResults : MetaMorpheusEngineResults
    {

        #region Private Fields

        private string v;

        #endregion Private Fields

        #region Public Constructors

        public MyErroredResults(MetaMorpheusEngine s, string v) : base(s)
        {
            this.v = v;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append(v);
            return sb.ToString();
        }

        #endregion Public Methods

    }
}