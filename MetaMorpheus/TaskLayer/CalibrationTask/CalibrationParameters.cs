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
            QValueCutoffForCalibratingPSMs = 0.005;
            WriteIndexedMzml = true;
        }

        public bool WriteIntermediateFiles { get; set; }
        public bool WriteIndexedMzml { get; set; }

        public int MinMS1IsotopicPeaksNeededForConfirmedIdentification { get; set; }
        public int MinMS2IsotopicPeaksNeededForConfirmedIdentification { get; set; }
        public int NumFragmentsNeededForEveryIdentification { get; set; }

        public double QValueCutoffForCalibratingPSMs { get; set; }
    }
}