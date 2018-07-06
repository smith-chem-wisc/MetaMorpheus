namespace EngineLayer.Calibration
{
    public class LabeledDataPoint
    {
        public readonly double ExperimentalMz;
        public readonly double RetentionTime;
        public readonly double LogTotalIonCurrent;
        public readonly double LogInjectionTime;
        public readonly double LogIntensity;
        public readonly double TheoreticalMz;
        public readonly double AbsoluteMzError;
        public readonly PeptideSpectralMatch Identification;

        public LabeledDataPoint(double experimentalMz, double rt, double logTotalIonCurrent, double logInjectionTime, double logIntensity, double theoreticalMz, PeptideSpectralMatch identification)
        {
            this.ExperimentalMz = experimentalMz;
            this.RetentionTime = rt;
            this.LogTotalIonCurrent = logTotalIonCurrent;
            this.LogInjectionTime = logInjectionTime;
            this.LogIntensity = logIntensity;
            this.TheoreticalMz = theoreticalMz;
            this.Identification = identification;
            this.AbsoluteMzError = experimentalMz - theoreticalMz;
        }
    }
}