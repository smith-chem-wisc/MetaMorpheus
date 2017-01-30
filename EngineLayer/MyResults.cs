using System;
using System.Text;

namespace EngineLayer
{
    public abstract class MyResults
    {

        #region Internal Fields

        internal TimeSpan Time;

        #endregion Internal Fields

        #region Protected Constructors

        protected MyResults(MyEngine s)
        {
            this.MyEngine = s;
        }

        #endregion Protected Constructors

        #region Protected Properties

        protected MyEngine MyEngine { get; private set; }

        protected abstract string StringForOutput { get; }

        #endregion Protected Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            if (MyEngine.Level <= 1)
            {
                sb.AppendLine(GetType().Name + ":");
                var ok = StringForOutput;
                if (!string.IsNullOrEmpty(ok))
                    sb.AppendLine(ok);
                sb.AppendLine("\tTime to run: " + Time);
            }
            else
            {
                sb.AppendLine("\t" + GetType().Name + ":");
                var ok = StringForOutput;
                if (!string.IsNullOrEmpty(ok))
                    sb.AppendLine(ok);
                sb.Append("\t\tTime to run: " + Time);
            }
            return sb.ToString();
        }

        #endregion Public Methods

    }
}