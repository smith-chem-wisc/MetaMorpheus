namespace EngineLayer.Util
{
    internal class ScanWithIndexAndNotchInfo(int notch, int scanIndex)
    {
        public readonly int Notch = notch;
        public readonly int ScanIndex = scanIndex;
    }

    /// <summary>
    /// For use in MiniClassicSearchEngine
    /// </summary>
    /// <param name="scan"></param>
    /// <param name="notch"></param>
    /// <param name="scanIndex"></param>
    internal class ExtendedScanWithIndexAndNotchInfo(Ms2ScanWithSpecificMass scan, int notch, int scanIndex) : ScanWithIndexAndNotchInfo(notch, scanIndex)
    {
        public readonly Ms2ScanWithSpecificMass TheScan = scan;
    }
}