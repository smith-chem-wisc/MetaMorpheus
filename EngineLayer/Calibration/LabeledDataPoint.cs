namespace EngineLayer.Calibration
{
    public class LabeledDataPoint
    {
        #region Public Fields

        public readonly double experimentalMz;
        public readonly double rt;
        public readonly double logTotalIonCurrent;
        public readonly double logInjectionTime;
        public readonly double logIntensity;
        public readonly double theoreticalMz;
        public readonly double absoluteMzError;
        public readonly PeptideSpectralMatch identification;

        #endregion Public Fields

        #region Protected Constructors

        public LabeledDataPoint(double experimentalMz, double rt, double logTotalIonCurrent, double logInjectionTime, double logIntensity, double theoreticalMz, PeptideSpectralMatch identification)
        {
            this.experimentalMz = experimentalMz;
            this.rt = rt;
            this.logTotalIonCurrent = logTotalIonCurrent;
            this.logInjectionTime = logInjectionTime;
            this.logIntensity = logIntensity;
            this.theoreticalMz = theoreticalMz;
            this.identification = identification;
            this.absoluteMzError = experimentalMz - theoreticalMz;
        }

        #endregion Protected Constructors
    }
}