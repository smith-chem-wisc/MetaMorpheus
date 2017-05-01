using Chemistry;
using MassSpectrometry;

namespace EngineLayer
{
    public class Ms2ScanWithSpecificMass
    {

        #region Public Constructors

        public Ms2ScanWithSpecificMass(IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> b, double precursorMz, int precursorCharge)
        {
            PrecursorMz = precursorMz;
            PrecursorCharge = precursorCharge;
            PrecursorMass = precursorMz.ToMass(precursorCharge);

            TheScan = b;

            OneBasedScanNumber = b.OneBasedScanNumber;
            OneBasedPrecursorScanNumber = b.OneBasedPrecursorScanNumber;

            RetentionTime = b.RetentionTime;

            NumPeaks = b.MassSpectrum.Size;

            TotalIonCurrent = b.TotalIonCurrent;
        }

        #endregion Public Constructors

        #region Public Properties

        public IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> TheScan { get; }
        public int OneBasedScanNumber { get; }
        public int OneBasedPrecursorScanNumber { get; }
        public double RetentionTime { get; }
        public int NumPeaks { get; }
        public double TotalIonCurrent { get; }
        public double PrecursorMz { get; private set; }
        public double PrecursorMass { get; private set; }
        public int PrecursorCharge { get; private set; }

        #endregion Public Properties

    }
}