using System;
using System.Text;

namespace EngineLayer
{
    public class MetaMorpheusEngineResults
    {
        #region Internal Fields

        internal TimeSpan Time;

        #endregion Internal Fields

        #region Public Constructors

        public MetaMorpheusEngineResults(MetaMorpheusEngine s)
        {
            MyEngine = s;
        }

        #endregion Public Constructors

        #region Protected Properties

        protected MetaMorpheusEngine MyEngine { get; }

        #endregion Protected Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(MyEngine.GetType().Name + " " + MyEngine.GetId());
            sb.Append("Time to run: " + Time);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}