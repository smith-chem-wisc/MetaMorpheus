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