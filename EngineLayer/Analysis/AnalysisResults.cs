using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Analysis
{
    public class AnalysisResults : MyResults
    {

        #region Private Fields

        private string output;

        #endregion Private Fields

        #region Public Constructors

        public AnalysisResults(AnalysisEngine s) : base(s)
        {
            output = "";
        }

        #endregion Public Constructors

        #region Public Properties

        public List<NewPsmWithFdr>[] AllResultingIdentifications { get; set; }
        public List<ProteinGroup>[] ProteinGroups { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(base.ToString());
            sb.AppendLine("All PSMS within 1% FDR: " + string.Join(", ", AllResultingIdentifications.Select(b => b.Count(c => c.qValue <= 0.01))));

            if (ProteinGroups != null && ProteinGroups.Any(s => s != null))
            {
                var numProteinsList = new List<int>();
                for (int i = 0; i < ProteinGroups.Length; i++)
                {
                    if (ProteinGroups[i] == null)
                        numProteinsList.Add(0);
                    else
                        numProteinsList.Add(ProteinGroups[i].Count(c => c.QValue <= 0.01));
                }

                sb.AppendLine("All proteins within 1% FDR: " + string.Join(", ", numProteinsList));
            }
            sb.Append(output);

            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal void AddText(string v)
        {
            output += v + Environment.NewLine;
        }

        #endregion Internal Methods

    }
}