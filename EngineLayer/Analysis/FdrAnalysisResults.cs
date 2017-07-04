using System.Collections.Generic;
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
            //sb.AppendLine("All PSMS within 1% FDR: " + string.Join(", ", AllResultingIdentifications.Select(b => b.Count(c => c.FdrInfo.QValue <= 0.01))));

            return sb.ToString();
        }

        #endregion Public Methods

    }
}