using MzLibUtil;

namespace EngineLayer
{
    public class CalibrationParameters
    {
        #region Public Constructors

        public CalibrationParameters()
        {
            PrecursorMassTolerance = new PpmTolerance(10);
            NonLinearCalibration = true;
            WriteIntermediateFiles = false;
            minMS1isotopicPeaksNeededForConfirmedIdentification = 3;
            minMS2isotopicPeaksNeededForConfirmedIdentification = 2;
            numFragmentsNeededForEveryIdentification = 6;
        }

        #endregion Public Constructors

        #region Public Properties

        public Tolerance PrecursorMassTolerance { get; set; }
        public bool NonLinearCalibration { get; set; }
        public bool WriteIntermediateFiles { get; set; }


        public int minMS1isotopicPeaksNeededForConfirmedIdentification { get; set; }
        public int minMS2isotopicPeaksNeededForConfirmedIdentification { get; set; }
        public int numFragmentsNeededForEveryIdentification { get; set; }

        #endregion Public Properties
    }
}