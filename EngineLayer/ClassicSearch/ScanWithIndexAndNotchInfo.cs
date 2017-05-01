namespace EngineLayer.ClassicSearch
{
    internal class ScanWithIndexAndNotchInfo
    {

        #region Public Fields

        public Ms2ScanWithSpecificMass theScan;
        public int notch;
        public int scanIndex;

        #endregion Public Fields

        #region Public Constructors

        public ScanWithIndexAndNotchInfo(Ms2ScanWithSpecificMass theScan, int notch, int scanIndex)
        {
            this.theScan = theScan;
            this.notch = notch;
            this.scanIndex = scanIndex;
        }

        #endregion Public Constructors

    }
}