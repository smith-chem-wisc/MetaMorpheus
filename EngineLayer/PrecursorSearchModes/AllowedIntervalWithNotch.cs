using MzLibUtil;

namespace EngineLayer
{
    public class AllowedIntervalWithNotch
    {
        #region Public Fields

        public DoubleRange allowedInterval;

        #endregion Public Fields

        #region Public Constructors

        public AllowedIntervalWithNotch(DoubleRange doubleRange, int j)
        {
            this.allowedInterval = doubleRange;
            this.Notch = j;
        }

        #endregion Public Constructors

        #region Public Properties

        public int Notch { get; }

        #endregion Public Properties
    }
}