using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class ProteinParsimonyResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public ProteinParsimonyResults(ProteinParsimonyEngine proteinAnalysisEngine) : base(proteinAnalysisEngine)
        {
        }

        #endregion Public Constructors

        #region Public Properties

        public List<ProteinGroup> ProteinGroups { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            if (ProteinGroups != null && ProteinGroups.Any(s => s != null))
            {
                var numProteinsList = new List<int>();
                if (ProteinGroups == null)
                    numProteinsList.Add(0);
                else
                {
                    int j = ProteinGroups.FindIndex(c => (c.QValue >= 0.01));
                    numProteinsList.Add(j);
                }

                sb.Append("All proteins within 1% FDR: " + string.Join(", ", numProteinsList));
            }
            return sb.ToString();
        }

        #endregion Public Methods
    }
}