using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Analysis
{
    public class AnalysisResults : MyResults
    {

        #region Public Constructors

        public AnalysisResults(AnalysisEngine s, List<NewPsmWithFdr>[] allResultingIdentifications, List<ProteinGroup> proteinGroups) : base(s)
        {
            this.AllResultingIdentifications = allResultingIdentifications;
            this.ProteinGroups = proteinGroups;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<NewPsmWithFdr>[] AllResultingIdentifications { get; private set; }
        public List<ProteinGroup> ProteinGroups { get; private set; }

        #endregion Public Properties

        #region Protected Properties

        protected override string StringForOutput
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append("\t\tAll PSMS within 1% FDR: " + string.Join(", ", AllResultingIdentifications.Select(b => b.Count(c => c.qValue <= 0.01))));

                if (ProteinGroups != null)
                sb.Append("\n\t\tAll proteins within 1% FDR: " + string.Join(", ", ProteinGroups.Count(c => ((c.QValue <= 0.01) && (c.isDecoy == false)))));

                return sb.ToString();
            }
        }

        #endregion Protected Properties

    }
}
