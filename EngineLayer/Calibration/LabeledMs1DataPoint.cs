namespace EngineLayer.Calibration
{
    public class LabeledMs1DataPoint : LabeledDataPoint
    {
        #region Public Constructors

        public LabeledMs1DataPoint(double mz, double rt, double LOGtotalIonCurrent, double LOGinjectionTime, double LOGintensity, double expectedMz, PeptideSpectralMatch identification)
            : base(mz, rt, LOGtotalIonCurrent, LOGinjectionTime, LOGintensity, expectedMz, identification)
        {
        }

        #endregion Public Constructors

        #region Public Properties

        public static string TabSeparatedHeader { get { return "mz\trt\tLOGTIC\tLOGInjectionTime\tLOGIntensity\tLabelTh\tlabelPPM"; } }

        #endregion Public Properties

        #region Public Methods

        public string Values()
        {
            return mz + "\t" + rt + "\t" + logTotalIonCurrent + "\t" + logInjectionTime + "\t" + logIntensity + "\t" + LabelTh + "\t" + LabelPpm;
        }

        #endregion Public Methods
    }
}