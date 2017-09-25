using MzLibUtil;
using System.Collections.Generic;


namespace EngineLayer
{
    public class OpenLowTheoSearchMode : MassDiffAcceptor
    {
        #region Public Constructors

        public OpenLowTheoSearchMode() : base("OpenLow")
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            if (scanPrecursorMass > peptideMass - 1)
                return 0;
            else
                return -1;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            yield return new AllowedIntervalWithNotch(new DoubleRange(-1, peptideMonoisotopicMass + 1), 0);
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

