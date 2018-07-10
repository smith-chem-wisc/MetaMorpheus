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
        }

        public bool WriteIntermediateFiles { get; set; }

        public int MinMS1IsotopicPeaksNeededForConfirmedIdentification { get; set; }
        public int MinMS2IsotopicPeaksNeededForConfirmedIdentification { get; set; }
        public int NumFragmentsNeededForEveryIdentification { get; set; }
    }
}