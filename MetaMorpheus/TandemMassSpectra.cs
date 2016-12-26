using System;
using System.Collections.Generic;

namespace MetaMorpheus
{
    public abstract class TandemMassSpectra : List<TandemMassSpectrum>
    {
        public string filename;

        public TandemMassSpectra() : base()
        {
        }

        public abstract void Load(string mzmlFilepath, int minimumAssumedPrecursorChargeState, int maximumAssumedPrecursorChargeState,
            int maximumNumberOfPeaks, bool assignChargeStates, bool deisotope, MassTolerance isotopicMZTolerance);

        public event EventHandler ReportTaskWithoutProgress;

        protected virtual void OnReportTaskWithoutProgress(EventArgs e)
        {
            ReportTaskWithoutProgress?.Invoke(null, e);
        }

        public event EventHandler ReportTaskWithProgress;

        protected virtual void OnReportTaskWithProgress(EventArgs e)
        {
            ReportTaskWithProgress?.Invoke(null, e);
        }

        public event EventHandler<int> UpdateProgress;

        protected virtual void OnUpdateProgress(int e)
        {
            UpdateProgress?.Invoke(null, e);
        }

        public IEnumerable<Tuple<TandemMassSpectrum, int>> GetTandemMassSpectraInMassRanges(double precursorMass, IEnumerable<double> acceptedPrecursorMassErrors, MassTolerance precursorMassTolerance)
        {
            foreach (double accepted_precursor_mass_error in acceptedPrecursorMassErrors)
            {
                foreach (Tuple<TandemMassSpectrum, int> spectrum in GetTandemMassSpectraInMassRange(precursorMass + accepted_precursor_mass_error, precursorMassTolerance))
                {
                    yield return spectrum;
                }
            }
        }

        public IEnumerable<Tuple<TandemMassSpectrum, int>> GetTandemMassSpectraInMassRange(double precursorMass, MassTolerance precursorMassTolerance)
        {
            double minimum_precursor_mass = precursorMass - precursorMassTolerance;
            double maximum_precursor_mass = precursorMass + precursorMassTolerance;

            int index = BinarySearch(NextHigherDouble(maximum_precursor_mass));
            if (index == Count)
            {
                index--;
            }
            while (index >= 0 && this[index].PrecursorMass >= minimum_precursor_mass)
            {
                if (this[index].PrecursorMass <= maximum_precursor_mass)
                {
                    yield return new Tuple<TandemMassSpectrum, int>(this[index], index);
                }
                index--;
            }
        }

        private int BinarySearch(double precursorMass)
        {
            int low_index = 0;
            int high_index = Count - 1;
            while (low_index <= high_index)
            {
                int mid_index = low_index + ((high_index - low_index) / 2);
                int comparison = this[mid_index].PrecursorMass.CompareTo(precursorMass);
                if (comparison == 0)
                {
                    return mid_index;
                }
                if (comparison < 0)
                {
                    low_index = mid_index + 1;
                }
                else
                {
                    high_index = mid_index - 1;
                }
            }
            return low_index;
        }

        // only works for positive doubles; does not handle special cases
        private static double NextHigherDouble(double value)
        {
            long bits = BitConverter.DoubleToInt64Bits(value);
            long next_double_bits = bits + 1;
            return BitConverter.Int64BitsToDouble(next_double_bits);
        }

        protected static List<MSPeak> FilterPeaks(List<MSPeak> peaks, int maximumNumberOfPeaks)
        {
            List<MSPeak> filtered_peaks = new List<MSPeak>(peaks);

            if (maximumNumberOfPeaks > 0 && filtered_peaks.Count > maximumNumberOfPeaks)
            {
                filtered_peaks.Sort(MSPeak.DescendingIntensityComparison);
                filtered_peaks.RemoveRange(maximumNumberOfPeaks, filtered_peaks.Count - maximumNumberOfPeaks);
            }

            return filtered_peaks;
        }
    }
}