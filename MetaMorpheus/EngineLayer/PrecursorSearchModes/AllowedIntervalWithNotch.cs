using MzLibUtil;

namespace EngineLayer
{
    public class AllowedIntervalWithNotch
    {
        public DoubleRange AllowedInterval;

        public AllowedIntervalWithNotch(DoubleRange doubleRange, int j)
        {
            AllowedInterval = doubleRange;
            Notch = j;
        }

        public int Notch { get; }
    }
}