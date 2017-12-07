using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TaskLayer
{
    public class MyTaskResults
    {
        #region Public Fields

        public List<string> newSpectra;//used else where, should be public
        public List<DbForTask> newDatabases;
        public TimeSpan Time;

        #endregion Public Fields

        #region Private Fields

        private readonly List<string> resultTexts;

        private readonly StringBuilder niceText = new StringBuilder();

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
            sb.AppendLine("Time to run task: " + Time);
            sb.AppendLine();
            sb.AppendLine();
            sb.AppendLine("--------------------------------------------------");
            if ((newSpectra != null && newSpectra.Any()) || (newDatabases != null && newDatabases.Any()))
            {
                sb.AppendLine();
                sb.AppendLine();
                sb.AppendLine("New files:");
                if (newSpectra != null && newSpectra.Any())
                {
                    sb.AppendLine("New spectra: ");
                    sb.AppendLine();
                    sb.AppendLine(string.Join(Environment.NewLine + "\t", newSpectra));
                }
                if (newDatabases != null && newDatabases.Any())
                {
                    sb.AppendLine("New databases: ");
                    sb.AppendLine(string.Join(Environment.NewLine + "\t", newDatabases.Select(b => b.FilePath)).ToString());
                }
                sb.AppendLine();
                sb.AppendLine();
                sb.AppendLine("--------------------------------------------------");
            }
            sb.AppendLine();
            sb.AppendLine();
            sb.AppendLine(niceText.ToString());
            sb.AppendLine();
            sb.AppendLine();
            sb.AppendLine("--------------------------------------------------");
            sb.AppendLine();
            sb.AppendLine();
            sb.AppendLine("Engine Results:");
            sb.AppendLine();
            foreach (var ok in resultTexts)
            {
                sb.AppendLine(ok);
                sb.AppendLine();
            }
            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal void AddResultText(string resultsText)
        {
            resultTexts.Add(resultsText);
        }

        internal void AddNiceText(string niceTextString)
        {
            niceText.AppendLine(niceTextString);
        }

        #endregion Internal Methods
    }
}