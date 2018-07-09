using MzLibUtil;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class RealTimeParameters
    {
        public RealTimeParameters()
        {
            DecoyType = DecoyType.Reverse;
            SearchTarget = true;
            MaxFragmentSize = 30000.0;
        }

        public DecoyType DecoyType { get; set; }
        public bool SearchTarget { get; set; }
        public double MaxFragmentSize { get; set; }
    }
}