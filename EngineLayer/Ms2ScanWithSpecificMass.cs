using Chemistry;
using MassSpectrometry;

namespace EngineLayer
{
    public class Ms2ScanWithSpecificMass
    {

        #region Public Constructors

        public Ms2ScanWithSpecificMass(IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> mzLibScan, IMzPeak precursorMonoisotopicPeak, int precursorCharge, string fullFilePath)
        {
            PrecursorMonoisotopicPeak = precursorMonoisotopicPeak;
            PrecursorCharge = precursorCharge;
            PrecursorMass = PrecursorMonoisotopicPeak.Mz.ToMass(precursorCharge);
            FullFilePath = fullFilePath;

            TheScan = mzLibScan;
        }

        #endregion Public Constructors

        #region Public Properties

        public IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> TheScan { get; }
        public IMzPeak PrecursorMonoisotopicPeak { get; }
        public double PrecursorMass { get; }
        public int PrecursorCharge { get; }
        public string FullFilePath { get; }

        #endregion Public Properties

    }
}