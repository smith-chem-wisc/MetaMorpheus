namespace TaskLayer
{
    public class CalibrationParameters
    {
        #region Public Constructors

        public CalibrationParameters()
        {
            WriteIntermediateFiles = false;
            MinMS1IsotopicPeaksNeededForConfirmedIdentification = 3;
            MinMS2IsotopicPeaksNeededForConfirmedIdentification = 2;
            NumFragmentsNeededForEveryIdentification = 10;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool WriteIntermediateFiles { get; set; }

        public int MinMS1IsotopicPeaksNeededForConfirmedIdentification { get; set; }
        public int MinMS2IsotopicPeaksNeededForConfirmedIdentification { get; set; }
        public int NumFragmentsNeededForEveryIdentification { get; set; }

        #endregion Public Properties
    }
}