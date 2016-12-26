using MSFileReaderLib;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Threading.Tasks;

namespace MetaMorpheus
{
    public class ThermoTandemMassSpectra : TandemMassSpectra
    {
        private const PrecursorMassType PRECURSOR_MASS_TYPE = PrecursorMassType.Monoisotopic;
        private const bool GET_PRECURSOR_MZ_AND_INTENSITY_FROM_MS1 = true;
        private const bool ALWAYS_USE_PRECURSOR_CHARGE_STATE_RANGE = false;

        public void Load(string rawFilepath, int minimumAssumedPrecursorChargeState, int maximumAssumedPrecursorChargeState,
            int maximumNumberOfPeaks,
            bool assignChargeStates)
        {
            OnReportTaskWithoutProgress(new EventArgs());

            IXRawfile2 raw = (IXRawfile2)new MSFileReader_XRawfile();

            raw.Open(rawFilepath);
            raw.SetCurrentController(0, 1);

            Dictionary<int, double[,]> ms1s;
            if (GET_PRECURSOR_MZ_AND_INTENSITY_FROM_MS1)
            {
                ms1s = new Dictionary<int, double[,]>();
            }

            int first_scan_number = -1;
            raw.GetFirstSpectrumNumber(ref first_scan_number);
            int last_scan_number = -1;
            raw.GetLastSpectrumNumber(ref last_scan_number);

            OnReportTaskWithProgress(new EventArgs());
            object progress_lock = new object();
            int spectra_processed = 0;
            int old_progress = 0;

            ParallelOptions parallel_options = new ParallelOptions();
            Parallel.For(first_scan_number, last_scan_number + 1, scan_number =>
            {
                string scan_filter = null;
                raw.GetFilterForScanNum(scan_number, ref scan_filter);

                if (!scan_filter.Contains(" ms "))
                {
                    string spectrum_id = "controllerType=0 controllerNumber=1 scan=" + scan_number.ToString();

                    double retention_time_minutes = double.NaN;
                    raw.RTFromScanNum(scan_number, ref retention_time_minutes);

                    int polarity = DeterminePolarity(scan_filter);

                    string fragmentation_method = GetFragmentationMethod(scan_filter);

                    double precursor_mz;
                    double precursor_intensity;
                    GetPrecursor(ms1s, raw, scan_number, scan_filter, first_scan_number, out precursor_mz, out precursor_intensity);

                    int charge = DeterminePrecursorCharge(raw, scan_number);
                    if (polarity < 0)
                    {
                        charge = -charge;
                    }

                    double[,] label_data = GetFragmentationData(raw, scan_number, scan_filter);
                    List<MSPeak> peaks = new List<MSPeak>(label_data.GetLength(1));
                    for (int peak_index = label_data.GetLowerBound(1); peak_index <= label_data.GetUpperBound(1); peak_index++)
                    {
                        peaks.Add(new MSPeak(label_data[(int)RawLabelDataColumn.MZ, peak_index],
                            label_data[(int)RawLabelDataColumn.Intensity, peak_index],
                            assignChargeStates && (int)RawLabelDataColumn.Charge < label_data.GetLength(0) ? (int)label_data[(int)RawLabelDataColumn.Charge, peak_index] : 0, polarity));
                    }

                    peaks = FilterPeaks(peaks, maximumNumberOfPeaks);
                    if (peaks.Count > 0)
                    {
                        for (int c = (ALWAYS_USE_PRECURSOR_CHARGE_STATE_RANGE || charge == 0 ? minimumAssumedPrecursorChargeState : charge);
                            c <= (ALWAYS_USE_PRECURSOR_CHARGE_STATE_RANGE || charge == 0 ? maximumAssumedPrecursorChargeState : charge); c++)
                        {
                            double precursor_mass = MSPeak.MassFromMZ(precursor_mz, c);

                            TandemMassSpectrum spectrum = new TandemMassSpectrum(rawFilepath, scan_number, spectrum_id, null, retention_time_minutes, fragmentation_method, precursor_mz, precursor_intensity, c, precursor_mass, peaks);

                            if (!SpectraLimits.limit || SpectraLimits.Contains(spectrum.SpectrumNumber))
                            {
                                lock (this)
                                {
                                    Add(spectrum);
                                }
                            }
                        }
                    }
                }

                lock (progress_lock)
                {
                    spectra_processed++;
                    int new_progress = (int)((double)spectra_processed / (last_scan_number - first_scan_number + 1) * 100);
                    if (new_progress > old_progress)
                    {
                        OnUpdateProgress(new_progress);
                        old_progress = new_progress;
                    }
                }
            });

            raw.Close();
        }

        private static string GetFragmentationMethod(string scanFilter)
        {
            if (scanFilter.Contains("cid"))
            {
                return "collision-induced dissociation";
            }
            else if (scanFilter.Contains("mpd"))
            {
                return "infrared multiphoton dissociation";
            }
            else if (scanFilter.Contains("pqd"))
            {
                return "pulsed q dissociation";
            }
            else if (scanFilter.Contains("hcd"))
            {
                return "high-energy collision-induced dissociation";
            }
            else if (scanFilter.Contains("ecd"))
            {
                return "electron capture dissociation";
            }
            else if (scanFilter.Contains("etd"))
            {
                return "electron transfer dissociation";
            }
            else
            {
                return null;
            }
        }

        private static void GetPrecursor(IDictionary<int, double[,]> ms1s, IXRawfile2 raw, int scanNumber, string scanFilter, int firstScanNumber, out double mz, out double intensity)
        {
            if (PRECURSOR_MASS_TYPE == PrecursorMassType.Monoisotopic)
            {
                mz = GetMonoisotopicMZ(raw, scanNumber, scanFilter);
            }
            if (GET_PRECURSOR_MZ_AND_INTENSITY_FROM_MS1)
            {
                GetAccurateMZAndIntensity(ms1s, raw, scanNumber, firstScanNumber, ref mz, out intensity);
            }
        }

        private static bool GetAccurateMZAndIntensity(IDictionary<int, double[,]> ms1s, IXRawfile2 raw, int scanNumber, int firstScanNumber, ref double mz, out double intensity)
        {
            scanNumber--;
            while (scanNumber >= firstScanNumber)
            {
                string scan_filter = null;
                raw.GetFilterForScanNum(scanNumber, ref scan_filter);

                if (scan_filter.Contains(" ms "))
                {
                    double[,] ms1;
                    lock (ms1s)
                    {
                        if (!ms1s.TryGetValue(scanNumber, out ms1))
                        {
                            if (scan_filter.Contains("FTMS"))
                            {
                                object labels_obj = null;
                                object flags_obj = null;
                                raw.GetLabelData(ref labels_obj, ref flags_obj, ref scanNumber);
                                ms1 = (double[,])labels_obj;
                            }
                            else
                            {
                                double centroid_peak_width = double.NaN;
                                object mass_list = null;
                                object peak_flags = null;
                                int array_size = -1;
                                raw.GetMassListFromScanNum(scanNumber, null, 0, 0, 0, 1, ref centroid_peak_width, ref mass_list, ref peak_flags, ref array_size);
                                ms1 = (double[,])mass_list;
                            }
                            ms1s.Add(scanNumber, ms1);
                        }
                    }

                    int index = -1;
                    for (int i = ms1.GetLowerBound(1); i <= ms1.GetUpperBound(1); i++)
                    {
                        if (index < 0 || Math.Abs(ms1[0, i] - mz) < Math.Abs(ms1[0, index] - mz))
                        {
                            index = i;
                        }
                    }
                    if (index >= 0)
                    {
                        mz = ms1[0, index];
                        intensity = ms1[1, index];
                        return true;
                    }
                    else
                    {
                        break;
                    }
                }
                else
                {
                    scanNumber--;
                }
            }

            intensity = double.NaN;
            return false;
        }

        private static double[,] GetFragmentationData(IXRawfile2 raw, int scanNumber, string scanFilter)
        {
            double[,] data;
            if (scanFilter.Contains("FTMS"))
            {
                object labels_obj = null;
                object flags_obj = null;
                raw.GetLabelData(ref labels_obj, ref flags_obj, ref scanNumber);
                data = (double[,])labels_obj;
            }
            else
            {
                double centroid_peak_width = double.NaN;
                object mass_list = null;
                object peak_flags = null;
                int array_size = -1;
                raw.GetMassListFromScanNum(scanNumber, null, 0, 0, 0, 1, ref centroid_peak_width, ref mass_list, ref peak_flags, ref array_size);
                data = (double[,])mass_list;
            }

            return data;
        }

        private static int DeterminePolarity(string scanFilter)
        {
            if (scanFilter.Contains(" + "))
            {
                return 1;
            }
            else if (scanFilter.Contains(" - "))
            {
                return -1;
            }
            else
            {
                throw new ArgumentException("Unknown polarity.");
            }
        }

        private static int DeterminePrecursorCharge(IXRawfile2 raw, int scanNumber)
        {
            object trailer_labels_obj = null;
            object trailer_values_obj = null;
            int trailer_array_size = -1;
            raw.GetTrailerExtraForScanNum(scanNumber, ref trailer_labels_obj, ref trailer_values_obj, ref trailer_array_size);
            string[] trailer_labels = (string[])trailer_labels_obj;
            string[] trailer_values = (string[])trailer_values_obj;

            int charge = -1;
            for (int trailer_index = trailer_labels.GetLowerBound(0); trailer_index <= trailer_labels.GetUpperBound(0); trailer_index++)
            {
                if (trailer_labels[trailer_index].StartsWith("Charge"))
                {
                    charge = int.Parse(trailer_values[trailer_index]);
                }
            }
            return charge;
        }

        private static double GetIsolationMZ(string scanFilter)
        {
            string temp_scan_filter = scanFilter.Substring(0, scanFilter.LastIndexOf('@'));
            double isolation_mz = double.Parse(temp_scan_filter.Substring(temp_scan_filter.LastIndexOf(' ') + 1), CultureInfo.InvariantCulture);

            return isolation_mz;
        }

        private static double GetMonoisotopicMZ(IXRawfile2 raw, int scanNumber, string scanFilter)
        {
            object labels_obj = null;
            object values_obj = null;
            int array_size = -1;
            raw.GetTrailerExtraForScanNum(scanNumber, ref labels_obj, ref values_obj, ref array_size);
            string[] labels = (string[])labels_obj;
            string[] values = (string[])values_obj;
            for (int i = labels.GetLowerBound(0); i <= labels.GetUpperBound(0); i++)
            {
                if (labels[i].StartsWith("Monoisotopic M/Z"))
                {
                    double monoisotopic_mz = double.Parse(values[i], CultureInfo.InvariantCulture);
                    if (monoisotopic_mz > 0.0)
                    {
                        return monoisotopic_mz;
                    }
                    else
                    {
                        break;
                    }
                }
            }

            return GetIsolationMZ(scanFilter);
        }

        public override void Load(string rawFilepath, int minimumAssumedPrecursorChargeState, int maximumAssumedPrecursorChargeState, int maximumNumberOfPeaks, bool assignChargeStates, bool deisotope, MassTolerance isotopicMZTolerance)
        {
            Load(rawFilepath, minimumAssumedPrecursorChargeState, maximumAssumedPrecursorChargeState,
                maximumNumberOfPeaks,
                assignChargeStates);
        }
    }

    internal enum PrecursorMassType
    {
        Isolation,
        Monoisotopic
    }

    internal enum RawLabelDataColumn
    {
        MZ = 0,
        Intensity = 1,
        Resolution = 2,
        NoiseBaseline = 3,
        NoiseLevel = 4,
        Charge = 5
    }
}