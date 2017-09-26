using MzLibUtil;

namespace TaskLayer
{
    public class CalibrationParameters
    {
        #region Public Constructors

        public CalibrationParameters()
        {
            PrecursorMassTolerance = new PpmTolerance(10);
            NonLinearCalibration = true;
            WriteIntermediateFiles = false;
            MinMS1isotopicPeaksNeededForConfirmedIdentification = 3;
            MinMS2isotopicPeaksNeededForConfirmedIdentification = 2;
            NumFragmentsNeededForEveryIdentification = 6;
        }

        #endregion Public Constructors

        #region Public Properties

        public Tolerance PrecursorMassTolerance { get; set; }
        public bool NonLinearCalibration { get; set; }
        public bool WriteIntermediateFiles { get; set; }


        public int MinMS1isotopicPeaksNeededForConfirmedIdentification { get; set; }
        public int MinMS2isotopicPeaksNeededForConfirmedIdentification { get; set; }
        public int NumFragmentsNeededForEveryIdentification { get; set; }

        #endregion Public Properties
    }
}