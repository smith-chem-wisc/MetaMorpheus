using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class ProteinAnalysisResults : MetaMorpheusEngineResults
    {

        #region Public Constructors

        public ProteinAnalysisResults(ProteinAnalysisEngine proteinAnalysisEngine) : base(proteinAnalysisEngine)
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