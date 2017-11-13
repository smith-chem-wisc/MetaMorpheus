using Chemistry;
using MassSpectrometry;

namespace EngineLayer
{
    public class Ms2ScanWithSpecificMass : IScan
    {
        #region Public Constructors

        public Ms2ScanWithSpecificMass(IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> mzLibScan, double precursorMonoisotopicPeakMz, int precursorCharge, string fullFilePath)
        {
            PrecursorMonoisotopicPeakMz = precursorMonoisotopicPeakMz;
            PrecursorCharge = precursorCharge;
            PrecursorMass = PrecursorMonoisotopicPeakMz.ToMass(precursorCharge);
            FullFilePath = fullFilePath;

            TheScan = mzLibScan;
        }

        #endregion Public Constructors

        #region Public Properties

        public IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> TheScan { get; }
        public double PrecursorMonoisotopicPeakMz { get; }
        public double PrecursorMass { get; }
        public int PrecursorCharge { get; }
        public string FullFilePath { get; }

        public int OneBasedScanNumber => TheScan.OneBasedScanNumber;

        public int? OneBasedPrecursorScanNumber => TheScan.OneBasedPrecursorScanNumber;

        public double RetentionTime => TheScan.RetentionTime;

        public int NumPeaks => TheScan.MassSpectrum.Size;

        public double TotalIonCurrent => TheScan.TotalIonCurrent;

        #endregion Public Properties
    }
}