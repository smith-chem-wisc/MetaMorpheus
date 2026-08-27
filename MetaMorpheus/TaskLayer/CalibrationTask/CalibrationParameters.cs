namespace TaskLayer
{
    public class CalibrationParameters
    {
        public CalibrationParameters()
        {
            WriteIntermediateFiles = false;
            MinMS1IsotopicPeaksNeededForConfirmedIdentification = 3;
            MinMS2IsotopicPeaksNeededForConfirmedIdentification = 2;
            NumFragmentsNeededForEveryIdentification = 10;
            QValueCutoffForCalibratingPSMs = 0.01;
            WriteIndexedMzml = true;
            OutputFormat = SpectraFileOutputFormat.MzML;
            SearchType = SearchType.Classic;
        }

        public bool WriteIntermediateFiles { get; set; }
        public bool WriteIndexedMzml { get; set; }

        /// <summary>
        /// Container for the calibrated spectra. mzML unless deliberately changed; see
        /// <see cref="SpectraFileOutputFormat"/> for the trade-off between the two MGF variants.
        /// WriteIndexedMzml is ignored for either of them.
        /// </summary>
        public SpectraFileOutputFormat OutputFormat { get; set; }

        public SearchType SearchType { get; set; }

        public int MinMS1IsotopicPeaksNeededForConfirmedIdentification { get; set; }
        public int MinMS2IsotopicPeaksNeededForConfirmedIdentification { get; set; }
        public int NumFragmentsNeededForEveryIdentification { get; set; }

        public double QValueCutoffForCalibratingPSMs { get; set; }
    }
}
