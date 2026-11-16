namespace EngineLayer.Util
{
    /// <summary>
    /// A scan in a candidate peptide's precursor window and the notch it fell in. A struct: the Classic search makes one
    /// for every peptide-scan pair, and each used to be an object for the garbage collector.
    /// </summary>
    internal readonly struct ScanWithIndexAndNotchInfo(int notch, int scanIndex)
    {
        public readonly int Notch = notch;
        public readonly int ScanIndex = scanIndex;
    }
}