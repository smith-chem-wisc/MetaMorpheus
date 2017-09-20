namespace EngineLayer.Calibration
{
    public class LabeledMs2DataPoint : LabeledDataPoint
    {
        #region Public Constructors

        public LabeledMs2DataPoint(double mz, double rt, double LOGtotalIonCurrent, double LOGinjectionTime, double LOGintensity, double label, Psm identification)
            : base(mz, rt, LOGtotalIonCurrent, LOGinjectionTime, LOGintensity, label, identification)
        {
        }

        #endregion Public Constructors

        #region Public Properties

        public static string TabSeparatedHeader { get { return "mz\trt\tLOGTIC\tLOGInjectionTime\tLOGIntensity\tLabel"; } }

        #endregion Public Properties

        #region Public Methods

        public string Values()
        {
            return mz + "\t" + rt + "\t" + logTotalIonCurrent + "\t" + logInjectionTime + "\t" + logIntensity + "\t" + LabelTh + "\t" + LabelPPM;
        }

        #endregion Public Methods
    }
}