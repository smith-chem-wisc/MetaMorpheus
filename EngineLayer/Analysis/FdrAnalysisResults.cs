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

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(base.ToString());
            sb.AppendLine("All target PSMS within 1% FDR: " + string.Join(", ", allResultingPeptides.Select(b => b.Count(c => c.FdrInfo.QValue <= 0.01 && !c.Pli.IsDecoy))));
            return sb.ToString();
        }

        #endregion Public Methods

    }
}