using MzLibUtil;
using System;
using System.Collections.Generic;


namespace EngineLayer
{
    public class OpenHighTheoSearchMode : MassDiffAcceptor
    {
        #region Public Constructors

        public OpenHighTheoSearchMode() : base("OpenHigh")
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            if (scanPrecursorMass < peptideMass + 1)
                return 0;
            else
                return -1;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            yield return new AllowedIntervalWithNotch(new DoubleRange(peptideMonoisotopicMass-1, Double.PositiveInfinity), 0);
        }

        public override string ToProseString()
        {
            return ("unboundedHigh");
        }

        public override string ToString()
        {
            return FileNameAddition + " OpenHighSearch";
        }

        #endregion Public Methods
    }
}

