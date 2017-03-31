using EngineLayer;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TaskLayer
{
    public class MyTaskResults
    {

        #region Public Fields

        public List<string> newSpectra;
        public List<DbForTask> newDatabases;
        public TimeSpan Time;

        #endregion Public Fields

        #region Private Fields

        private readonly List<string> resultTexts;

        #endregion Private Fields

        #region Internal Constructors

        internal MyTaskResults(MetaMorpheusTask s)
        {
            resultTexts = new List<string>();
        }

        #endregion Internal Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(GetType().Name);
            sb.AppendLine("Time to run: " + Time);
            if (newSpectra != null && newSpectra.Any())
                sb.AppendLine("New spectra: " + string.Join(", ", newSpectra));
            if (newDatabases != null && newDatabases.Any())
                sb.AppendLine("New databases: " + string.Join(", ", newDatabases.Select(b => b.FileName)));
            foreach (var ok in resultTexts)
                sb.Append(ok);
            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal void AddResultText(MetaMorpheusEngineResults indexResults)
        {
            resultTexts.Add(indexResults.ToString());
        }

        #endregion Internal Methods

    }
}