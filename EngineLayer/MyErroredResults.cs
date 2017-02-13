using System.Text;

namespace EngineLayer
{
    public class MyErroredResults : MyResults
    {

        #region Private Fields

        private string v;

        #endregion Private Fields

        #region Public Constructors

        public MyErroredResults(MyEngine s, string v) : base(s)
        {
            this.v = v;
        }

        #endregion Public Constructors

        #region Protected Properties


        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.Append(v);
            return sb.ToString();
        }

        #endregion Protected Properties

    }
}