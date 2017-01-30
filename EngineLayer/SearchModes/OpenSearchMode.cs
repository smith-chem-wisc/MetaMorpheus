using Spectra;
using System.Collections.Generic;

namespace EngineLayer
{
    public class OpenSearchMode : SearchMode
    {

        #region Public Constructors

        public OpenSearchMode() : base("OpenSearch")
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

        #endregion Public Methods

    }
}