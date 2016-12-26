using System.Collections.Generic;

namespace MetaMorpheus
{
    public class TandemMassSpectrum
    {
        public string Filename { get; private set; }

        public int SpectrumNumber { get; private set; }

        public string SpectrumId { get; private set; }

        public string SpectrumTitle { get; private set; }

        public double RetentionTimeMinutes { get; private set; }

        public string FragmentationMethod { get; private set; }

        public double PrecursorMZ { get; private set; }

        public double PrecursorIntensity { get; private set; }

        public int PrecursorCharge { get; private set; }

        public double PrecursorMass { get; private set; }

        public List<double> mzs { get; private set; }

        public double[] Intensities { get; private set; }

        public double TotalIntensity { get; private set; }

        public TandemMassSpectrum(string filename, int spectrumNumber, string spectrumId, string spectrumTitle, double retentionTimeMinutes, string fragmentationMethod, double precursorMZ, double precursorIntensity, int precursorCharge, double precursorMass, List<MSPeak> peaks)
        {
            Filename = filename;
            SpectrumNumber = spectrumNumber;
            SpectrumId = spectrumId;
            SpectrumTitle = spectrumTitle;
            RetentionTimeMinutes = retentionTimeMinutes;
            FragmentationMethod = fragmentationMethod;
            PrecursorMZ = precursorMZ;
            PrecursorIntensity = precursorIntensity;
            PrecursorCharge = precursorCharge;
            PrecursorMass = precursorMass;
            TotalIntensity = 0.0;
            if (peaks != null)
            {
                peaks.Sort(MSPeak.AscendingMassComparison);
                mzs = new List<double>();
                Intensities = new double[peaks.Count];
                for (int p = 0; p < peaks.Count; p++)
                {
                    MSPeak peak = peaks[p];
                    mzs.Add(peak.mz);
                    Intensities[p] = peak.Intensity;
                    TotalIntensity += peak.Intensity;
                }
            }
        }

        public static int AscendingPrecursorMassComparison(TandemMassSpectrum left, TandemMassSpectrum right)
        {
            return left.PrecursorMass.CompareTo(right.PrecursorMass);
        }
    }
}