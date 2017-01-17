using Spectra;
using System.Collections.Generic;

namespace InternalLogicEngineLayer
{
    public class OpenSearchMode : SearchMode
    {

        #region Public Constructors

        public OpenSearchMode(string s) : base(s)
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override bool Accepts(double scanPrecursorMass, double peptideMass)
        {
            return true;
        }

        public override IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            yield return new DoubleRange(double.MinValue, double.MaxValue);
        }

        public override string SearchModeString()
        {
            return "OpenSearch";
        }

        #endregion Public Methods

    }
}