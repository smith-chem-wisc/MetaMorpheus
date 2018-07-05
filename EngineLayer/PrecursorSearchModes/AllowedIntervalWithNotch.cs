using MzLibUtil;

namespace EngineLayer
{
    public class AllowedIntervalWithNotch
    {
        public DoubleRange allowedInterval;

        public AllowedIntervalWithNotch(DoubleRange doubleRange, int j)
        {
            this.allowedInterval = doubleRange;
            this.Notch = j;
        }

        public int Notch { get; }
    }
}