namespace EngineLayer.Calibration
{
    public abstract class LabeledDataPoint
    {
        #region Public Fields

        public readonly double mz;
        public readonly double rt;
        public readonly double logTotalIonCurrent;
        public readonly double logInjectionTime;

        public readonly double logIntensity;

        public readonly double expectedMZ;

        public readonly Psm identification;

        #endregion Public Fields

        #region Public Constructors

        protected LabeledDataPoint(double mz, double rt, double logTotalIonCurrent, double logInjectionTime, double logIntensity, double expectedMZ, Psm identification)
        {
            this.mz = mz;
            this.rt = rt;
            this.logTotalIonCurrent = logTotalIonCurrent;
            this.logInjectionTime = logInjectionTime;

            this.logIntensity = logIntensity;

            this.expectedMZ = expectedMZ;

            this.identification = identification;
        }

        #endregion Public Constructors

        #region Public Properties

        public double LabelTh { get { return mz - expectedMZ; } }
        public double LabelPpm { get { return (mz - expectedMZ) / (expectedMZ) * 1e6; } }

        #endregion Public Properties
    }
}