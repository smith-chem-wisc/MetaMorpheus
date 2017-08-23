using MzLibUtil;

namespace EngineLayer
{
    public class CalibrationParameters
    {
        public Tolerance PrecursorMassTolerance { get; set; }
        public bool NonLinearCalibration { get; set; }
        public bool WriteIntermediateFiles { get; set; }
    }
}
