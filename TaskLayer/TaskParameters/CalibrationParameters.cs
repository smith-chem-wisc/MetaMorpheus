using MzLibUtil;

namespace EngineLayer
{
    public class CalibrationParameters
    {
        public CalibrationParameters()
        {
            PrecursorMassTolerance = new PpmTolerance(10);
            NonLinearCalibration = true;
            WriteIntermediateFiles = false;
        }
        public Tolerance PrecursorMassTolerance { get; set; }
        public bool NonLinearCalibration { get; set; }
        public bool WriteIntermediateFiles { get; set; }
    }
}
