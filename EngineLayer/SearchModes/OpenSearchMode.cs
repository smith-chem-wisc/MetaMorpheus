using MzLibUtil;
using System;
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

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            return 0;
        }

        public override IEnumerable<Tuple<DoubleRange, int>> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            yield return new Tuple<DoubleRange, int>(new DoubleRange(double.MinValue, double.MaxValue), 0);
        }

        public override string ToString()
        {
            return FileNameAddition + " OpenSearch";
        }

        #endregion Public Methods

    }
}