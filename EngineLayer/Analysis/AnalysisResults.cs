using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Analysis
{
    public class AnalysisResults : MyResults
    {

        #region Public Constructors

        public AnalysisResults(AnalysisEngine s, List<NewPsmWithFdr>[] allResultingIdentifications, List<ProteinGroup>[] proteinGroups) : base(s)
        {
            this.AllResultingIdentifications = allResultingIdentifications;
            this.ProteinGroups = proteinGroups;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<NewPsmWithFdr>[] AllResultingIdentifications { get; private set; }
        public List<ProteinGroup>[] ProteinGroups { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("All PSMS within 1% FDR: " + string.Join(", ", AllResultingIdentifications.Select(b => b.Count(c => c.qValue <= 0.01))));

            var check = ProteinGroups.Where(s => s != null);
            if (check.Any())
            {
                var numProteinsList = new List<int>();
                for (int i = 0; i < ProteinGroups.Length; i++)
                {
                    if (ProteinGroups[i] == null)
                        numProteinsList.Add(0);
                    else
                        numProteinsList.Add(ProteinGroups[i].Count(c => c.QValue <= 0.01));
                }

                sb.Append("All proteins within 1% FDR: " + string.Join(", ", numProteinsList));
            }

            return sb.ToString();
        }

        #endregion Public Methods

    }
}