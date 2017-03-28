using System.Text;

namespace TaskLayer
{
    public class EverythingRunnerResults
    {

        #region Private Fields

        private readonly string status;

        #endregion Private Fields

        #region Public Constructors

        public EverythingRunnerResults(string status)
        {
            this.status = status;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(status);
            return sb.ToString();
        }

        #endregion Public Methods

    }
}