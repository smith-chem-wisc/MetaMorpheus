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
            MyEngine = s;
        }

        #endregion Protected Constructors

        #region Protected Properties

        protected MyEngine MyEngine { get; }

        #endregion Protected Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(GetType().Name);
            sb.AppendLine("Time to run: " + Time);
            return sb.ToString();
        }

        #endregion Public Methods

    }
}