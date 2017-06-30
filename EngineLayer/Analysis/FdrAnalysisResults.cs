using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Analysis
{
    public class FdrAnalysisResults : MetaMorpheusEngineResults
    {
        #region Public Fields

        public List<PsmParent>[] allResultingPeptides;

        #endregion Public Fields

        #region Internal Fields

        internal Dictionary<string, int>[] allModsOnPeptides;
        internal Dictionary<string, int>[] allModsSeen;

        #endregion Internal Fields

        #region Public Constructors

        public FdrAnalysisResults(FdrAnalysisEngine s) : base(s)
        {
        }

        #endregion Public Constructors

        #region Public Properties

        public List<ProteinGroup>[] ProteinGroups { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            //sb.AppendLine("All PSMS within 1% FDR: " + string.Join(", ", AllResultingIdentifications.Select(b => b.Count(c => c.FdrInfo.QValue <= 0.01))));

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
                sb.AppendLine(string.Join(Environment.NewLine, allModsSeen[i].OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
                sb.AppendLine("Search mode " + i + " Mods in database:");
                sb.AppendLine(string.Join(Environment.NewLine, allModsOnPeptides[i].OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            }

            return sb.ToString();
        }

        #endregion Public Methods
    }
}