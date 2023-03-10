namespace EngineLayer
{
    internal class ScanWithIndexAndNotchInfo
    {
        public Ms2ScanWithSpecificMass TheScan;
        public int Notch;
        public int ScanIndex;

        public ScanWithIndexAndNotchInfo(Ms2ScanWithSpecificMass theScan, int notch, int scanIndex)
        {
            TheScan = theScan;
            Notch = notch;
            ScanIndex = scanIndex;
        }
    }
}