using InternalLogicEngineLayer;
using System.Text;

namespace InternalLogicTaskLayer
{
    internal class EverythingRunnerResults : MyResults
    {

        #region Public Constructors

        public EverythingRunnerResults(MyEngine s) : base(s)
        {
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override string GetStringForOutput()
        {
            var sb = new StringBuilder();
            return sb.ToString();
        }

        #endregion Protected Methods

    }
}