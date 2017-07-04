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
            sb.Append(base.ToString());
            //sb.AppendLine("All PSMS within 1% FDR: " + string.Join(", ", AllResultingIdentifications.Select(b => b.Count(c => c.FdrInfo.QValue <= 0.01))));

            if (ProteinGroups != null && ProteinGroups.Any(s => s != null))
            {
                sb.AppendLine();
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

                sb.Append("All proteins within 1% FDR: " + string.Join(", ", numProteinsList));
            }

            return sb.ToString();
        }

        #endregion Public Methods

    }
}