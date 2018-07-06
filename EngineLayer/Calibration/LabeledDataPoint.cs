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

        public LabeledDataPoint(double experimentalMz, double rt, double logTotalIonCurrent, double logInjectionTime,
            double logIntensity, double theoreticalMz, PeptideSpectralMatch identification)
        {
            ExperimentalMz = experimentalMz;
            RetentionTime = rt;
            LogTotalIonCurrent = logTotalIonCurrent;
            LogInjectionTime = logInjectionTime;
            LogIntensity = logIntensity;
            TheoreticalMz = theoreticalMz;
            Identification = identification;
            AbsoluteMzError = experimentalMz - theoreticalMz;
        }
    }
}