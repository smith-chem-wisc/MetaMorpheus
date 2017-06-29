using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class MetaMorpheusEngineResults
    {

        #region Internal Fields

        internal TimeSpan Time;

        #endregion Internal Fields

        #region Private Fields

        private List<string> outputStrings = new List<string>();

        #endregion Private Fields

        #region Protected Constructors

        public MetaMorpheusEngineResults(MetaMorpheusEngine s)
        {
            MyEngine = s;
        }

        #endregion Protected Constructors

        #region Protected Properties

        protected MetaMorpheusEngine MyEngine { get; }

        #endregion Protected Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(GetType().Name);
            sb.AppendLine("Time to run: " + Time);
            sb.Append(string.Join(Environment.NewLine, outputStrings));
            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal void AddText(string v)
        {
            outputStrings.Add(v);
        }

        #endregion Internal Methods

    }
}