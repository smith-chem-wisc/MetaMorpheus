using MassSpectrometry;

namespace EngineLayer
{
    public interface IScan
    {
        #region Public Properties

        string FullFilePath { get; }
        int OneBasedScanNumber { get; }
        int OneBasedPrecursorScanNumber { get; }
        double RetentionTime { get; }
        int NumPeaks { get; }
        double TotalIonCurrent { get; }
        int PrecursorCharge { get; }
        IMzPeak PrecursorMonoisotopicPeak { get; }
        double PrecursorMass { get; }

        #endregion Public Properties
    }
}