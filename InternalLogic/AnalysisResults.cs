using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class AnalysisResults : MyResults
    {
        #region Public Constructors

        public AnalysisResults(AnalysisEngine s, List<NewPsmWithFDR>[] allResultingIdentifications, List<ProteinGroup> proteinGroups) : base(s)
        {
            this.allResultingIdentifications = allResultingIdentifications;
            this.proteinGroups = proteinGroups;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<NewPsmWithFDR>[] allResultingIdentifications { get; private set; }
        public List<ProteinGroup> proteinGroups { get; private set; }

        #endregion Public Properties

        #region Protected Methods

        protected override string GetStringForOutput()
        {
            var sb = new StringBuilder();
            sb.Append("\t\tAll PSMS within 1% FDR: " + string.Join(", ", allResultingIdentifications.Select(b => b.Count(c => c.QValue <= 0.01))));

            if (proteinGroups != null)
                sb.Append("\t\tAll proteins within 1% FDR: " + string.Join(", ", proteinGroups.Select(proteinGroup => proteinGroup.QValue <= 0.01)));

            return sb.ToString();
        }

        #endregion Protected Methods
    }
}