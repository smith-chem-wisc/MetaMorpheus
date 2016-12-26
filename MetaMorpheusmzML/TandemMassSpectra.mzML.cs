using Ionic.Zlib;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;
using System.Xml;
using System.Xml.XPath;

namespace MetaMorpheus
{
    public class MzMLTandemMassSpectra : TandemMassSpectra
    {
        private const bool GET_PRECURSOR_MZ_AND_INTENSITY_FROM_MS1 = true;
        private const bool ALWAYS_USE_PRECURSOR_CHARGE_STATE_RANGE = false;

        internal static void ReadDataFromSpectrumNavigator(XPathNodeIterator binaryDataArrayChildNodes, out double[] mz, out double[] intensity)
        {
            mz = null;
            intensity = null;

            int word_length_in_bytes = 0;
            bool zlib_compressed = false;
            ArrayDataType array_data_type = ArrayDataType.Unknown;
            foreach (XPathNavigator navigator in binaryDataArrayChildNodes)
            {
                if (navigator.Name.Equals("cvParam", StringComparison.OrdinalIgnoreCase))
                {
                    if (navigator.GetAttribute("name", string.Empty).Equals("32-bit float", StringComparison.OrdinalIgnoreCase))
                    {
                        word_length_in_bytes = 4;
                    }
                    else if (navigator.GetAttribute("name", string.Empty).Equals("64-bit float", StringComparison.OrdinalIgnoreCase))
                    {
                        word_length_in_bytes = 8;
                    }
                    else if (navigator.GetAttribute("name", string.Empty).Equals("zlib compression", StringComparison.OrdinalIgnoreCase))
                    {
                        zlib_compressed = true;
                    }
                    else if (navigator.GetAttribute("name", string.Empty).Equals("m/z array", StringComparison.OrdinalIgnoreCase))
                    {
                        array_data_type = ArrayDataType.MZ;
                    }
                    else if (navigator.GetAttribute("name", string.Empty).Equals("intensity array", StringComparison.OrdinalIgnoreCase))
                    {
                        array_data_type = ArrayDataType.Intensity;
                    }
                }
                else if (navigator.Name.Equals("binary", StringComparison.OrdinalIgnoreCase))
                {
                    if (array_data_type == ArrayDataType.MZ)
                    {
                        mz = ReadBase64EncodedDoubleArray(navigator.InnerXml, word_length_in_bytes, zlib_compressed);
                    }
                    else if (array_data_type == ArrayDataType.Intensity)
                    {
                        intensity = ReadBase64EncodedDoubleArray(navigator.InnerXml, word_length_in_bytes, zlib_compressed);
                    }
                }
            }
        }

        private static double[] ReadBase64EncodedDoubleArray(string base64EncodedData, int wordLengthInBytes, bool zlibCompressed)
        {
            byte[] bytes = Convert.FromBase64String(base64EncodedData);
            if (zlibCompressed)
            {
                bytes = ZlibStream.UncompressBuffer(bytes);
            }
            double[] doubles = new double[bytes.Length / wordLengthInBytes];
            if (wordLengthInBytes == 4)
            {
                for (int i = doubles.GetLowerBound(0); i <= doubles.GetUpperBound(0); i++)
                {
                    doubles[i] = BitConverter.ToSingle(bytes, i * wordLengthInBytes);
                }
            }
            else if (wordLengthInBytes == 8)
            {
                Buffer.BlockCopy(bytes, 0, doubles, 0, bytes.Length);
            }
            return doubles;
        }

        private static List<MSPeak> AssignChargeStates(IList<MSPeak> peaks, int maxAbsoluteCharge, int polarity, MassTolerance isotopicMZTolerance)
        {
            List<MSPeak> new_peaks = new List<MSPeak>();

            for (int i = 0; i < peaks.Count; i++)
            {
                int j = i + 1;
                List<int> charges = new List<int>();
                while (j < peaks.Count)
                {
                    if (peaks[j].MZ > (peaks[i].MZ + Constants.C12_C13_MASS_DIFFERENCE) + isotopicMZTolerance)
                    {
                        break;
                    }

                    for (int c = polarity * maxAbsoluteCharge; polarity > 0 ? c >= 1 : c <= -1; c -= polarity)
                    {
                        if (Math.Abs(MassTolerance.CalculateMassError(peaks[j].MZ, peaks[i].MZ + Constants.C12_C13_MASS_DIFFERENCE / c, isotopicMZTolerance.Units)) <= isotopicMZTolerance.Value)
                        {
                            new_peaks.Add(new MSPeak(peaks[i].MZ, peaks[i].Intensity, c, polarity));
                            charges.Add(c);
                        }
                    }

                    j++;
                }
                if (charges.Count == 0)
                {
                    new_peaks.Add(new MSPeak(peaks[i].MZ, peaks[i].Intensity, 0, polarity));
                }
            }

            return new_peaks;
        }

        private static List<MSPeak> Deisotope(IEnumerable<MSPeak> peaks, int maxAbsoluteCharge, int polarity, MassTolerance isotopicMZTolerance)
        {
            List<MSPeak> new_peaks = new List<MSPeak>(peaks);

            int p = new_peaks.Count - 1;
            while (p >= 1)
            {
                int q = p - 1;
                bool removed = false;
                while (q >= 0)
                {
                    if (new_peaks[p].MZ > (new_peaks[q].MZ + Constants.C12_C13_MASS_DIFFERENCE) + isotopicMZTolerance)
                    {
                        break;
                    }

                    if (new_peaks[p].Intensity < new_peaks[q].Intensity)
                    {
                        if (polarity == 0)
                        {
                            if (Math.Abs(MassTolerance.CalculateMassError(new_peaks[p].MZ, new_peaks[q].MZ + Constants.C12_C13_MASS_DIFFERENCE, isotopicMZTolerance.Units)) <= isotopicMZTolerance.Value)
                            {
                                new_peaks.RemoveAt(p);
                                removed = true;
                                break;
                            }
                        }
                        else
                        {
                            for (int c = polarity; polarity > 0 ? c <= maxAbsoluteCharge : c >= -maxAbsoluteCharge; c += polarity)
                            {
                                if (Math.Abs(MassTolerance.CalculateMassError(new_peaks[p].MZ, new_peaks[q].MZ + Constants.C12_C13_MASS_DIFFERENCE / Math.Abs(c), isotopicMZTolerance.Units)) <= isotopicMZTolerance.Value)
                                {
                                    new_peaks.RemoveAt(p);
                                    removed = true;
                                    break;
                                }
                            }
                        }
                        if (removed)
                        {
                            break;
                        }
                    }

                    q--;
                }

                p--;
            }

            return new_peaks;
        }

        public override void Load(string mzmlFilepath, int minimumAssumedPrecursorChargeState, int maximumAssumedPrecursorChargeState, int maximumNumberOfPeaks, bool assignChargeStates, bool deisotope, MassTolerance isotopicMZTolerance)
        {
            OnReportTaskWithoutProgress(EventArgs.Empty);

            XmlDocument mzML_temp = new XmlDocument();
            mzML_temp.Load(mzmlFilepath);
            XPathNavigator mzML = mzML_temp.CreateNavigator();

            XmlNamespaceManager xnm = new XmlNamespaceManager(mzML.NameTable);
            xnm.AddNamespace("mzML", mzML_temp.DocumentElement.NamespaceURI);

            Dictionary<string, XPathNodeIterator> referenceable_param_groups = new Dictionary<string, XPathNodeIterator>();
            foreach (XPathNavigator referenceable_param_group in mzML.Select("//mzML:mzML/mzML:referenceableParamGroupList/mzML:referenceableParamGroup", xnm))
            {
                referenceable_param_groups.Add(referenceable_param_group.GetAttribute("id", string.Empty), referenceable_param_group.SelectChildren(XPathNodeType.All));
            }

            ParallelOptions parallel_options = new ParallelOptions();

            Dictionary<string, SpectrumNavigator> ms1s = null;
            if (GET_PRECURSOR_MZ_AND_INTENSITY_FROM_MS1)
            {
                ms1s = new Dictionary<string, SpectrumNavigator>();
                Parallel.ForEach(mzML.Select("//mzML:mzML/mzML:run/mzML:spectrumList/mzML:spectrum", xnm).Cast<XPathNavigator>(), parallel_options, spectrum_navigator =>

                {
                    string scan_id = spectrum_navigator.GetAttribute("id", string.Empty);
                    int ms_level = -1;

                    foreach (XPathNavigator spectrum_child_navigator in spectrum_navigator.SelectChildren(XPathNodeType.All))
                    {
                        if (spectrum_child_navigator.Name.Equals("cvParam", StringComparison.OrdinalIgnoreCase))
                        {
                            if (spectrum_child_navigator.GetAttribute("name", string.Empty).Equals("ms level", StringComparison.OrdinalIgnoreCase))
                            {
                                ms_level = int.Parse(spectrum_child_navigator.GetAttribute("value", string.Empty));
                            }
                        }
                        else if (spectrum_child_navigator.Name.Equals("referenceableParamGroupRef", StringComparison.OrdinalIgnoreCase))
                        {
                            foreach (XPathNavigator navigator in referenceable_param_groups[spectrum_child_navigator.GetAttribute("ref", string.Empty)])
                            {
                                if (navigator.Name.Equals("cvParam", StringComparison.OrdinalIgnoreCase))
                                {
                                    if (navigator.GetAttribute("name", string.Empty).Equals("ms level", StringComparison.OrdinalIgnoreCase))
                                    {
                                        ms_level = int.Parse(navigator.GetAttribute("value", string.Empty));
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    if (ms_level == 1)
                    {
                        lock (ms1s)
                        {
                            ms1s.Add(scan_id, new SpectrumNavigator(spectrum_navigator, xnm));
                        }
                    }
                }
                );
            }

            int num_spectra = int.Parse(mzML.SelectSingleNode("//mzML:mzML/mzML:run/mzML:spectrumList", xnm).GetAttribute("count", string.Empty));

            OnReportTaskWithProgress(EventArgs.Empty);
            object progress_lock = new object();
            int spectra_processed = 0;
            int old_progress = 0;

            Parallel.ForEach(mzML.Select("//mzML:mzML/mzML:run/mzML:spectrumList/mzML:spectrum", xnm).Cast<XPathNavigator>(), parallel_options, spectrum_navigator =>

            {
                int spectrum_index = int.Parse(spectrum_navigator.GetAttribute("index", string.Empty));
                int spectrum_number = spectrum_index + 1;
                string spectrum_id = spectrum_navigator.GetAttribute("id", string.Empty);
                string spectrum_title = null;
                int ms_level = -1;
                int polarity = 0;
                double retention_time_minutes = double.NaN;
                string precursor_scan_id = null;
                double precursor_mz = double.NaN;
                int charge = 0;
                double precursor_intensity = double.NaN;
                string fragmentation_method = "collision-induced dissociation";
                double[] mz = null;
                double[] intensity = null;

                foreach (XPathNavigator spectrum_child_navigator in spectrum_navigator.SelectChildren(XPathNodeType.All))
                {
                    if (spectrum_child_navigator.Name.Equals("cvParam", StringComparison.OrdinalIgnoreCase))
                    {
                        if (spectrum_child_navigator.GetAttribute("name", string.Empty).Equals("ms level", StringComparison.OrdinalIgnoreCase))
                        {
                            ms_level = int.Parse(spectrum_child_navigator.GetAttribute("value", string.Empty));
                        }
                        else if (spectrum_child_navigator.GetAttribute("name", string.Empty).Equals("positive scan", StringComparison.OrdinalIgnoreCase))
                        {
                            polarity = 1;
                        }
                        else if (spectrum_child_navigator.GetAttribute("name", string.Empty).Equals("negative scan", StringComparison.OrdinalIgnoreCase))
                        {
                            polarity = -1;
                        }
                        else if (spectrum_child_navigator.GetAttribute("name", string.Empty).Equals("spectrum title", StringComparison.OrdinalIgnoreCase))
                        {
                            spectrum_title = spectrum_child_navigator.GetAttribute("value", string.Empty);
                        }
                    }
                    else if (spectrum_child_navigator.Name.Equals("referenceableParamGroupRef", StringComparison.OrdinalIgnoreCase))
                    {
                        foreach (XPathNavigator navigator in referenceable_param_groups[spectrum_child_navigator.GetAttribute("ref", string.Empty)])
                        {
                            if (navigator.Name.Equals("cvParam", StringComparison.OrdinalIgnoreCase))
                            {
                                if (navigator.GetAttribute("name", string.Empty).Equals("ms level", StringComparison.OrdinalIgnoreCase))
                                {
                                    ms_level = int.Parse(navigator.GetAttribute("value", string.Empty));
                                }
                                else if (navigator.GetAttribute("name", string.Empty).Equals("positive scan", StringComparison.OrdinalIgnoreCase))
                                {
                                    polarity = 1;
                                }
                                else if (navigator.GetAttribute("name", string.Empty).Equals("negative scan", StringComparison.OrdinalIgnoreCase))
                                {
                                    polarity = -1;
                                }
                            }
                        }
                    }
                    else if (spectrum_child_navigator.Name.Equals("scanList", StringComparison.OrdinalIgnoreCase))
                    {
                        foreach (XPathNavigator navigator in spectrum_child_navigator.Select("mzML:scan/mzML:cvParam", xnm))
                        {
                            if (navigator.GetAttribute("name", string.Empty).Equals("scan start time", StringComparison.OrdinalIgnoreCase))
                            {
                                retention_time_minutes = double.Parse(navigator.GetAttribute("value", string.Empty), CultureInfo.InvariantCulture);
                                if (navigator.GetAttribute("unitName", string.Empty).StartsWith("s", StringComparison.OrdinalIgnoreCase))
                                {
                                    retention_time_minutes = TimeSpan.FromSeconds(retention_time_minutes).TotalMinutes;
                                }
                            }
                        }
                    }
                    else if (spectrum_child_navigator.Name.Equals("precursorList", StringComparison.OrdinalIgnoreCase))
                    {
                        XPathNavigator precursor_node = spectrum_child_navigator.SelectSingleNode("mzML:precursor", xnm);
                        precursor_scan_id = precursor_node.GetAttribute("spectrumRef", string.Empty);
                        foreach (XPathNavigator navigator in precursor_node.Select("mzML:selectedIonList/mzML:selectedIon/mzML:cvParam", xnm))
                        {
                            if (navigator.GetAttribute("name", string.Empty).Equals("selected ion m/z", StringComparison.OrdinalIgnoreCase))
                            {
                                precursor_mz = double.Parse(navigator.GetAttribute("value", string.Empty), CultureInfo.InvariantCulture);
                            }
                            else if (navigator.GetAttribute("name", string.Empty).Equals("charge state", StringComparison.OrdinalIgnoreCase))
                            {
                                charge = int.Parse(navigator.GetAttribute("value", string.Empty));
                                if (polarity < 0)
                                {
                                    charge = -charge;
                                }
                            }
                            else if (navigator.GetAttribute("name", string.Empty).Equals("peak intensity", StringComparison.OrdinalIgnoreCase))
                            {
                                precursor_intensity = double.Parse(navigator.GetAttribute("value", string.Empty), CultureInfo.InvariantCulture);
                            }
                        }
                        XPathNavigator navigator2 = spectrum_child_navigator.SelectSingleNode("mzML:precursor/mzML:activation/mzML:cvParam", xnm);
                        if (navigator2 != null)
                        {
                            fragmentation_method = navigator2.GetAttribute("name", string.Empty);
                        }
                    }
                    else if (spectrum_child_navigator.Name.Equals("binaryDataArrayList", StringComparison.OrdinalIgnoreCase))
                    {
                        ReadDataFromSpectrumNavigator(spectrum_child_navigator.Select("mzML:binaryDataArray/*", xnm), out mz, out intensity);
                    }
                    if (ms_level == 1)
                    {
                        break;
                    }
                }

                if (ms_level >= 2)
                {
                    if (GET_PRECURSOR_MZ_AND_INTENSITY_FROM_MS1 && precursor_scan_id != null)
                    {
                        SpectrumNavigator ms1;
                        if (ms1s.TryGetValue(precursor_scan_id, out ms1))
                        {
                            double[] ms1_mz;
                            double[] ms1_intensity;
                            ms1.GetSpectrum(out ms1_mz, out ms1_intensity);
                            int index = -1;
                            for (int i = ms1_mz.GetLowerBound(0); i <= ms1_mz.GetUpperBound(0); i++)
                            {
                                if (index < 0 || Math.Abs(ms1_mz[i] - precursor_mz) < Math.Abs(ms1_mz[index] - precursor_mz))
                                {
                                    index = i;
                                }
                            }
                            precursor_mz = ms1_mz[index];
                            precursor_intensity = ms1_intensity[index];
                        }
                    }

                    if (mz != null && intensity != null && mz.Length > 0 && intensity.Length > 0)
                    {
                        List<MSPeak> peaks = new List<MSPeak>(mz.Length);
                        for (int i = 0; i < mz.Length; i++)
                        {
                            peaks.Add(new MSPeak(mz[i], intensity[i], 0, polarity));
                        }

                        peaks = FilterPeaks(peaks, maximumNumberOfPeaks);
                        if (peaks.Count > 0)
                        {
                            peaks.Sort(MSPeak.AscendingMZComparison);
                            for (int c = (ALWAYS_USE_PRECURSOR_CHARGE_STATE_RANGE || charge == 0 ? minimumAssumedPrecursorChargeState : charge);
                                c <= (ALWAYS_USE_PRECURSOR_CHARGE_STATE_RANGE || charge == 0 ? maximumAssumedPrecursorChargeState : charge); c++)
                            {
                                List<MSPeak> new_peaks = peaks;
                                if (assignChargeStates)
                                {
                                    new_peaks = AssignChargeStates(new_peaks, c, polarity, isotopicMZTolerance);
                                    if (deisotope)
                                    {
                                        new_peaks = Deisotope(new_peaks, c, polarity, isotopicMZTolerance);
                                    }
                                }

                                double precursor_mass = MSPeak.MassFromMZ(precursor_mz, c);

                                TandemMassSpectrum spectrum = new TandemMassSpectrum(mzmlFilepath, spectrum_number, spectrum_id, spectrum_title, retention_time_minutes, fragmentation_method, precursor_mz, precursor_intensity, c, precursor_mass, new_peaks);
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
                }

                lock (progress_lock)
                {
                    spectra_processed++;
                    int new_progress = (int)((double)spectra_processed / num_spectra * 100);
                    if (new_progress > old_progress)
                    {
                        OnUpdateProgress(new_progress);
                        old_progress = new_progress;
                    }
                }
            }
        );
        }
    }

    internal enum ArrayDataType
    {
        Unknown,
        MZ,
        Intensity
    }

    internal class SpectrumNavigator
    {
        private XPathNavigator navigator;
        private XmlNamespaceManager xmlNamespaceManager;
        private double[] mz;
        private double[] intensity;

        internal SpectrumNavigator(XPathNavigator navigator, XmlNamespaceManager xmlNamespaceManager)
        {
            this.navigator = navigator;
            this.xmlNamespaceManager = xmlNamespaceManager;
        }

        internal void GetSpectrum(out double[] mz, out double[] intensity)
        {
            if (this.mz == null || this.intensity == null)
            {
                MzMLTandemMassSpectra.ReadDataFromSpectrumNavigator(navigator.Select("mzML:binaryDataArrayList/mzML:binaryDataArray/*", xmlNamespaceManager), out mz, out intensity);
                this.mz = mz;
                this.intensity = intensity;
            }
            else
            {
                mz = this.mz;
                intensity = this.intensity;
            }
        }
    }
}