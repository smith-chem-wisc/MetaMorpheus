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

        public LocalMS2Scan(IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> b, int charge)
        {
            TheScan = b;
            MonoisotopicPrecursorCharge = charge;
            MonoisotopicPrecursorMZ = TheScan.SelectedIonGuessMonoisotopicMZ;
            PrecursorMass = MonoisotopicPrecursorMZ.ToMass(MonoisotopicPrecursorCharge);

            OneBasedScanNumber = b.OneBasedScanNumber;

            RetentionTime = b.RetentionTime;

            NumPeaks = b.MassSpectrum.Size;

            double monoisotopicPrecursorIntensityhere = TheScan.SelectedIonGuessMonoisotopicIntensity;
            MonoisotopicPrecursorIntensity = monoisotopicPrecursorIntensityhere;

            TotalIonCurrent = b.TotalIonCurrent;
        }

        public LocalMS2Scan(double precursorMass)
        {
            this.PrecursorMass = precursorMass;
        }

        #endregion Public Constructors

        #region Public Properties

        public double PrecursorMass { get; private set; }

        #endregion Public Properties

        #region Internal Properties

        internal IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> TheScan { get; private set; }
        internal int OneBasedScanNumber { get; private set; }
        internal double RetentionTime { get; private set; }
        internal int MonoisotopicPrecursorCharge { get; private set; }
        internal int NumPeaks { get; private set; }
        internal double MonoisotopicPrecursorIntensity { get; private set; }
        internal double TotalIonCurrent { get; private set; }

        #endregion Internal Properties

    }
}