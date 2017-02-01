using Chemistry;
using MassSpectrometry;
using Spectra;
using System;

namespace EngineLayer
{
    public class LocalMS2Scan
    {

        #region Public Constructors

        public LocalMS2Scan(IMsDataScan<IMzSpectrum<MzPeak>> b, int charge)
        {
            TheScan = b;
            b.TryGetSelectedIonGuessMonoisotopicMZ(out MonoisotopicPrecursorMZ);
            MonoisotopicPrecursorCharge = charge;
            PrecursorMass = MonoisotopicPrecursorMZ.ToMass(MonoisotopicPrecursorCharge);

            OneBasedScanNumber = b.OneBasedScanNumber;

            RetentionTime = b.RetentionTime;

            NumPeaks = b.MassSpectrum.Count;

            double monoisotopicPrecursorIntensityhere;
            b.TryGetSelectedIonGuessMonoisotopicIntensity(out monoisotopicPrecursorIntensityhere);
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

        internal IMsDataScan<IMzSpectrum<MzPeak>> TheScan { get; private set; }
        internal readonly double MonoisotopicPrecursorMZ;
        internal int OneBasedScanNumber { get; private set; }
        internal double RetentionTime { get; private set; }
        internal int MonoisotopicPrecursorCharge { get; private set; }
        internal int NumPeaks { get; private set; }
        internal double MonoisotopicPrecursorIntensity { get; private set; }
        internal double TotalIonCurrent { get; private set; }

        #endregion Internal Properties

    }
}