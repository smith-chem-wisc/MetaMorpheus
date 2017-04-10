using Chemistry;
using MassSpectrometry;

namespace EngineLayer
{
    public class LocalMS2Scan
    {

        #region Internal Fields

        internal readonly double MonoisotopicPrecursorMZ;

        #endregion Internal Fields

        #region Public Constructors

        public LocalMS2Scan(IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> b)
        {
            TheScan = b;
            PrecursorCharge = b.SelectedIonGuessChargeStateGuess.Value;
            MonoisotopicPrecursorMZ = TheScan.SelectedIonGuessMonoisotopicMZ.Value;
            MonoisotopicPrecursorMass = MonoisotopicPrecursorMZ.ToMass(PrecursorCharge);

            OneBasedScanNumber = b.OneBasedScanNumber;
            OneBasedPrecursorScanNumber = b.OneBasedPrecursorScanNumber;

            RetentionTime = b.RetentionTime;

            NumPeaks = b.MassSpectrum.Size;

            MonoisotopicPrecursorIntensity = TheScan.SelectedIonGuessMonoisotopicIntensity.Value;

            TotalIonCurrent = b.TotalIonCurrent;
        }

        public LocalMS2Scan(double precursorMass)
        {
            this.MonoisotopicPrecursorMass = precursorMass;
        }

        #endregion Public Constructors

        #region Public Properties

        public double MonoisotopicPrecursorMass { get; private set; }

        #endregion Public Properties

        #region Internal Properties

        internal IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> TheScan { get; private set; }
        internal int OneBasedScanNumber { get; private set; }
        internal int OneBasedPrecursorScanNumber { get; private set; }
        internal double RetentionTime { get; private set; }
        internal int PrecursorCharge { get; private set; }
        internal int NumPeaks { get; private set; }
        internal double MonoisotopicPrecursorIntensity { get; private set; }
        internal double TotalIonCurrent { get; private set; }

        #endregion Internal Properties

    }
}