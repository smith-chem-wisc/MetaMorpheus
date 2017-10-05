namespace EngineLayer
{
    public class CalibrationParameters
    {
        #region Public Constructors

        public CalibrationParameters()
        {
            NonLinearCalibration = true;
            WriteIntermediateFiles = false;
        }

        #endregion Public Constructors

        #region Public Properties

        public bool NonLinearCalibration { get; set; }
        public bool WriteIntermediateFiles { get; set; }

        #endregion Public Properties
    }
}