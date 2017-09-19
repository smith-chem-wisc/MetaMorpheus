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

        public readonly double label;

        public readonly Psm identification;

        #endregion Public Fields

        #region Public Constructors

        public LabeledDataPoint(double mz, double rt, double logTotalIonCurrent, double logInjectionTime, double logIntensity, double label, Psm identification)
        {
            this.mz = mz;
            this.rt = rt;
            this.logTotalIonCurrent = logTotalIonCurrent;
            this.logInjectionTime = logInjectionTime;

            this.logIntensity = logIntensity;

            this.label = label;

            this.identification = identification;
        }

        #endregion Public Constructors
    }
}