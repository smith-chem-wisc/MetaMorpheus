using Chemistry;
using MassSpectrometry;

namespace EngineLayer
{
    public class Ms2ScanWithSpecificMass
    {

        #region Public Constructors

        public Ms2ScanWithSpecificMass(IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> b, double precursorMz, int precursorCharge, string fileNameWithoutExtension)
        {
            PrecursorMz = precursorMz;
            PrecursorCharge = precursorCharge;
            PrecursorMass = precursorMz.ToMass(precursorCharge);
            FileNameWithoutExtension = fileNameWithoutExtension;

            TheScan = b;
        }

        #endregion Public Constructors

        #region Public Properties

        public IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> TheScan { get; }
        public double PrecursorMz { get; }
        public double PrecursorMass { get; }
        public int PrecursorCharge { get; }
        public string FileNameWithoutExtension { get; }

        #endregion Public Properties

    }
}