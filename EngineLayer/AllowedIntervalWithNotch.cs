using MzLibUtil;

namespace EngineLayer
{
    public class AllowedIntervalWithNotch
    {

        #region Public Fields

        public DoubleRange allowedInterval;
        public int notch;

        #endregion Public Fields

        #region Public Constructors

        public AllowedIntervalWithNotch(DoubleRange doubleRange, int j)
        {
            this.allowedInterval = doubleRange;
            this.notch = j;
        }

        #endregion Public Constructors

    }
}