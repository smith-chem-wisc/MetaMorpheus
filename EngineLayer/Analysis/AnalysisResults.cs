using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Analysis
{
    public class AnalysisResults : MetaMorpheusEngineResults
    {

        #region Internal Fields

        internal Dictionary<string, int>[] allModsOnPeptides;
        internal Dictionary<string, int>[] allModsSeen;

        #endregion Internal Fields

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
            sb.AppendLine("All PSMS within 1% FDR: " + string.Join(", ", AllResultingIdentifications.Select(b => b.Count(c => c.QValue <= 0.01))));

            if (ProteinGroups != null && ProteinGroups.Any(s => s != null))
            {
                var numProteinsList = new List<int>();
                for (int i = 0; i < ProteinGroups.Length; i++)
                {
                    if (ProteinGroups[i] == null)
                        numProteinsList.Add(0);
                    else
                    {
                        int j = ProteinGroups[i].FindIndex(c => (c.QValue >= 0.01));
                        numProteinsList.Add(j);
                    }
                }

                sb.AppendLine("All proteins within 1% FDR: " + string.Join(", ", numProteinsList));
            }

            for (int i = 0; i < allModsOnPeptides.Length; i++)
            {
                sb.AppendLine("Search mode " + i + " Mods seen:");
                sb.AppendLine(string.Join(Environment.NewLine, allModsSeen[i].OrderBy(b => -b.Value).Select(b => b.Key + "\t" + b.Value)));
                sb.AppendLine("Search mode " + i + " Mods on proteins:");
                sb.AppendLine(string.Join(Environment.NewLine, allModsOnPeptides[i].OrderBy(b => -b.Value).Select(b => b.Key + "\t" + b.Value)));
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