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

        #region Public Properties

        public MetaMorpheusEngine MyEngine { get; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append("Time to run: " + Time);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}