using EngineLayer;
using System.Text;

namespace TaskLayer
{
    public class EverythingRunnerResults : MyResults
    {

        #region Public Constructors

        public EverythingRunnerResults(MyEngine s) : base(s)
        {
        }

        #endregion Public Constructors

        #region Protected Properties

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine( base.ToString());
            return sb.ToString();
        }

        #endregion Protected Properties

    }
}