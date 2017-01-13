using System;
using System.Text;

namespace InternalLogicEngineLayer
{
    public abstract class MyResults
    {
        #region Public Constructors

        public MyResults(MyEngine s)
        {
            this.s = s;
        }

        #endregion Public Constructors

        #region Public Properties

        internal TimeSpan Time;

        #endregion Public Properties

        #region Protected Properties

        protected MyEngine s { get; private set; }

        #endregion Protected Properties

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            if (s.Level <= 1)
            {
                sb.AppendLine(GetType().Name.ToString() + ":");
                var ok = GetStringForOutput();
                if (!string.IsNullOrEmpty(ok))
                    sb.AppendLine(ok);
                sb.AppendLine("\tTime to run: " + Time);
            }
            else
            {
                sb.AppendLine("\t" + GetType().Name.ToString() + ":");
                var ok = GetStringForOutput();
                if (!string.IsNullOrEmpty(ok))
                    sb.AppendLine(ok);
                sb.Append("\t\tTime to run: " + Time);
            }
            return sb.ToString();
        }

        #endregion Public Methods

        #region Protected Methods

        protected abstract string GetStringForOutput();

        #endregion Protected Methods
    }
}