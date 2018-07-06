namespace EngineLayer
{
    public interface IScan
    {
        string FullFilePath { get; }
        int OneBasedScanNumber { get; }
        int? OneBasedPrecursorScanNumber { get; }
        double RetentionTime { get; }
        int NumPeaks { get; }
        double TotalIonCurrent { get; }
        int PrecursorCharge { get; }
        double PrecursorMonoisotopicPeakMz { get; }
        double PrecursorMass { get; }
    }
}