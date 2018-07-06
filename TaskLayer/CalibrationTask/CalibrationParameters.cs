namespace TaskLayer
{
    public class CalibrationParameters
    {
        public bool WriteIntermediateFiles { get; set; }
        public int MinMS1IsotopicPeaksNeededForConfirmedIdentification { get; set; } = 3;
        public int MinMS2IsotopicPeaksNeededForConfirmedIdentification { get; set; } = 2;
        public int NumFragmentsNeededForEveryIdentification { get; set; } = 10;
    }
}