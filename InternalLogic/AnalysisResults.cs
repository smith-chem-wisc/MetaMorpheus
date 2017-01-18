using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InternalLogicEngineLayer
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

        #region Protected Methods

        protected override string StringForOutput
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append("\t\tAll PSMS within 1% FDR: " + string.Join(", ", AllResultingIdentifications.Select(b => b.Count(c => c.qValue <= 0.01))));

            if (proteinGroups != null)
                sb.Append("\n\t\tAll proteins within 1% FDR: " + string.Join(", ", proteinGroups.Count(c => c.QValue <= 0.01)));

                return sb.ToString();
            }
        }

        #endregion Protected Methods
    }
}