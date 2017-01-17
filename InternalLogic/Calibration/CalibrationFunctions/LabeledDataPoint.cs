namespace InternalLogicCalibration
{
    public class LabeledDataPoint
    {

        #region Public Fields

        public double[] inputs;
        public double output;

        #endregion Public Fields

        #region Public Constructors

        public LabeledDataPoint(double[] v1, double v2)
        {
            inputs = v1;
            output = v2;
        }

        #endregion Public Constructors

    }
}