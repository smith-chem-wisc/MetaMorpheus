namespace EngineLayer
{
    internal class ScanWithIndexAndNotchInfo
    {
        public Ms2ScanWithSpecificMass theScan;
        public int notch;
        public int scanIndex;

        public ScanWithIndexAndNotchInfo(Ms2ScanWithSpecificMass theScan, int notch, int scanIndex)
        {
            this.theScan = theScan;
            this.notch = notch;
            this.scanIndex = scanIndex;
        }
    }
}