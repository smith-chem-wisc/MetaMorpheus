using Chemistry;
using MassSpectrometry;
using Spectra;

namespace InternalLogicEngineLayer
{
    public class LocalMS2Scan
    {

        #region Public Constructors

        public LocalMS2Scan(IMsDataScan<IMzSpectrum<MzPeak>> b)
        {
            theScan = b;
            double monoisotopicPrecursorMZhere;
            b.TryGetSelectedIonGuessMonoisotopicMZ(out monoisotopicPrecursorMZhere);
            MonoisotopicPrecursorMZ = monoisotopicPrecursorMZhere;

            int monoisotopicPrecursorChargehere;
            b.TryGetSelectedIonGuessChargeStateGuess(out monoisotopicPrecursorChargehere);
            monoisotopicPrecursorCharge = monoisotopicPrecursorChargehere;

            PrecursorMass = MonoisotopicPrecursorMZ.ToMass(monoisotopicPrecursorCharge);

            OneBasedScanNumber = b.OneBasedScanNumber;

            RetentionTime = b.RetentionTime;

            NumPeaks = b.MassSpectrum.Count;

            double monoisotopicPrecursorIntensityhere;
            b.TryGetSelectedIonGuessMonoisotopicIntensity(out monoisotopicPrecursorIntensityhere);
            monoisotopicPrecursorIntensity = monoisotopicPrecursorIntensityhere;

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

        internal IMsDataScan<IMzSpectrum<MzPeak>> theScan { get; private set; }
        internal double MonoisotopicPrecursorMZ { get; private set; }
        internal int OneBasedScanNumber { get; private set; }
        internal double RetentionTime { get; private set; }
        internal int monoisotopicPrecursorCharge { get; private set; }
        internal int NumPeaks { get; private set; }
        internal double monoisotopicPrecursorIntensity { get; private set; }
        internal double TotalIonCurrent { get; private set; }

        #endregion Internal Properties

    }
}