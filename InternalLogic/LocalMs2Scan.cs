using Chemistry;
using MassSpectrometry;
using Spectra;
using System;

namespace InternalLogicEngineLayer
{
    public class LocalMs2Scan : IComparable, IComparable<double>
    {
        #region Public Constructors

        public LocalMs2Scan(IMsDataScan<IMzSpectrum<MzPeak>> b)
        {
            theScan = b;
            double monoisotopicPrecursorMZhere;
            b.TryGetSelectedIonGuessMonoisotopicMZ(out monoisotopicPrecursorMZhere);
            monoisotopicPrecursorMZ = monoisotopicPrecursorMZhere;

            int monoisotopicPrecursorChargehere;
            b.TryGetSelectedIonGuessChargeStateGuess(out monoisotopicPrecursorChargehere);
            monoisotopicPrecursorCharge = monoisotopicPrecursorChargehere;

            precursorMass = monoisotopicPrecursorMZ.ToMass(monoisotopicPrecursorCharge);

            OneBasedScanNumber = b.OneBasedScanNumber;

            RetentionTime = b.RetentionTime;

            NumPeaks = b.MassSpectrum.Count;

            double monoisotopicPrecursorIntensityhere;
            b.TryGetSelectedIonGuessMonoisotopicIntensity(out monoisotopicPrecursorIntensityhere);
            monoisotopicPrecursorIntensity = monoisotopicPrecursorIntensityhere;

            TotalIonCurrent = b.TotalIonCurrent;
        }

        #endregion Public Constructors

        #region Public Properties

        public double precursorMass { get; private set; }

        #endregion Public Properties

        #region Internal Properties

        internal IMsDataScan<IMzSpectrum<MzPeak>> theScan { get; private set; }
        internal double monoisotopicPrecursorMZ { get; private set; }
        internal int OneBasedScanNumber { get; private set; }
        internal double RetentionTime { get; private set; }
        internal int monoisotopicPrecursorCharge { get; private set; }
        internal int NumPeaks { get; private set; }
        internal double monoisotopicPrecursorIntensity { get; private set; }
        internal double TotalIonCurrent { get; private set; }

        #endregion Internal Properties

        #region Public Methods

        public int CompareTo(double other)
        {
            return precursorMass.CompareTo(other);
        }

        public int CompareTo(object obj)
        {
            var other = obj as LocalMs2Scan;
            if (other != null)
                return precursorMass.CompareTo(other.precursorMass);
            return precursorMass.CompareTo((double)obj);
        }

        #endregion Public Methods
    }
}