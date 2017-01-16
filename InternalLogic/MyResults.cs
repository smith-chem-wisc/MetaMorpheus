using System;
using System.Text;

namespace InternalLogicEngineLayer
{
    public abstract class MyResults
    {

        #region Internal Fields

        internal TimeSpan Time;

        #endregion Internal Fields

        #region Protected Constructors

        protected MyResults(MyEngine s)
        {
            this.s = s;
        }

        #endregion Protected Constructors

        #region Protected Properties

        protected MyEngine s { get; private set; }

        #endregion Protected Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append("MetaMorpheus version: " + MyEngine.MetaMorpheusVersion);
            if (s.Level <= 1)
            {
                sb.AppendLine(GetType().Name + ":");
                var ok = GetStringForOutput();
                if (!string.IsNullOrEmpty(ok))
                    sb.AppendLine(ok);
                sb.AppendLine("\tTime to run: " + Time);
            }
            else
            {
                sb.AppendLine("\t" + GetType().Name + ":");
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