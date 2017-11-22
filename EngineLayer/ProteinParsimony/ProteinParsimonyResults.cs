using System.Collections.Generic;
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
            return sb.ToString();
        }

        #endregion Public Methods
    }
}