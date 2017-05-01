using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Globalization;

namespace EngineLayer
{
    public class DaltonsAroundZeroWithInnerSearchMode : SearchMode
    {

        #region Private Fields

        private readonly SearchMode dotSearchMode;
        private readonly double daltonsAroundZero;

        #endregion Private Fields

        #region Public Constructors

        public DaltonsAroundZeroWithInnerSearchMode(SearchMode innerSearchMode, double daltonsAroundZero) : base("coisolation" + daltonsAroundZero.ToString("F3", CultureInfo.InvariantCulture) + innerSearchMode.FileNameAddition)
        {
            this.dotSearchMode = innerSearchMode;
            this.daltonsAroundZero = daltonsAroundZero;
            NumNotches = 1 + innerSearchMode.NumNotches;
        }

        #endregion Public Constructors

        #region Public Methods

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            if (Math.Abs(scanPrecursorMass - peptideMass) <= daltonsAroundZero)
                return 0;
            int innerResult = dotSearchMode.Accepts(scanPrecursorMass, peptideMass);
            return innerResult < 0 ? innerResult : innerResult + 1;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            foreach (var heh in dotSearchMode.GetAllowedPrecursorMassIntervals(peptideMonoisotopicMass))
            {
                yield return new AllowedIntervalWithNotch(heh.allowedInterval, heh.notch + 1);
            }
            yield return new AllowedIntervalWithNotch(new DoubleRange(peptideMonoisotopicMass - daltonsAroundZero, peptideMonoisotopicMass + daltonsAroundZero), 0);
        }

        #endregion Public Methods

    }
}