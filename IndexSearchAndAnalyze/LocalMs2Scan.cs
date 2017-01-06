using Chemistry;
using MassSpectrometry;
using Spectra;
using System;

namespace IndexSearchAndAnalyze
{
    internal class LocalMs2Scan : IComparable, IComparable<double>
    {
        public IMsDataScan<IMzSpectrum<MzPeak>> b { get; private set; }

        public LocalMs2Scan(IMsDataScan<IMzSpectrum<MzPeak>> b)
        {
            this.b = b;
            double monoisotopicPrecursorMZ;
            b.TryGetSelectedIonGuessMonoisotopicMZ(out monoisotopicPrecursorMZ);
            this.monoisotopicPrecursorMZ = monoisotopicPrecursorMZ;

            int monoisotopicPrecursorCharge;
            b.TryGetSelectedIonGuessChargeStateGuess(out monoisotopicPrecursorCharge);
            this.monoisotopicPrecursorCharge = monoisotopicPrecursorCharge;

            precursorMass = monoisotopicPrecursorMZ.ToMass(monoisotopicPrecursorCharge);

            OneBasedScanNumber = b.OneBasedScanNumber;

            RetentionTime = b.RetentionTime;

            NumPeaks = b.MassSpectrum.Count;

            double monoisotopicPrecursorIntensity;
            b.TryGetSelectedIonGuessMonoisotopicIntensity(out monoisotopicPrecursorIntensity);
            this.monoisotopicPrecursorIntensity = monoisotopicPrecursorIntensity;

            TotalIonCurrent = b.TotalIonCurrent;
        }

        public double precursorMass { get; private set; }
        public double monoisotopicPrecursorMZ { get; internal set; }
        public int OneBasedScanNumber { get; internal set; }
        public double RetentionTime { get; internal set; }
        public int monoisotopicPrecursorCharge { get; internal set; }
        public int NumPeaks { get; internal set; }
        public double monoisotopicPrecursorIntensity { get; internal set; }
        public double TotalIonCurrent { get; internal set; }

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
    }
}