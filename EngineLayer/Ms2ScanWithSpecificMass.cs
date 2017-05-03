using Chemistry;
using MassSpectrometry;

namespace EngineLayer
{
    public class Ms2ScanWithSpecificMass
    {

        #region Public Constructors

        public Ms2ScanWithSpecificMass(IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> b, IMzPeak precursorMonoisotopicPeak, int precursorCharge, string fileNameWithoutExtension)
        {
            PrecursorMonoisotopicPeak = precursorMonoisotopicPeak;
            PrecursorCharge = precursorCharge;
            PrecursorMass = PrecursorMonoisotopicPeak.Mz.ToMass(precursorCharge);
            FileNameWithoutExtension = fileNameWithoutExtension;

            TheScan = b;
        }

        #endregion Public Constructors

        #region Public Properties

        public IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> TheScan { get; }
        public IMzPeak PrecursorMonoisotopicPeak { get; }
        public double PrecursorMass { get; }
        public int PrecursorCharge { get; }
        public string FileNameWithoutExtension { get; }

        #endregion Public Properties

    }
}