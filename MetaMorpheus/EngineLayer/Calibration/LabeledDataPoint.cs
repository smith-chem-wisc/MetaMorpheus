﻿namespace EngineLayer.Calibration
{
    public class LabeledDataPoint
    {
        public readonly double ExperimentalMz;
        public readonly int ScanNumber;
        public readonly double LogTotalIonCurrent;
        public readonly double LogInjectionTime;
        public readonly double LogIntensity;
        public readonly double TheoreticalMz;
        public readonly double RelativeMzError;
        public readonly SpectralMatch Identification;

        public LabeledDataPoint(double experimentalMz, int scanNumber, double logTotalIonCurrent, double logInjectionTime, double logIntensity, double theoreticalMz, SpectralMatch identification)
        {
            this.ExperimentalMz = experimentalMz;
            this.ScanNumber = scanNumber;
            this.LogTotalIonCurrent = logTotalIonCurrent;
            this.LogInjectionTime = logInjectionTime;
            this.LogIntensity = logIntensity > 0 ? logIntensity : 0; //peaks with zero intensity can cause downstream crash
            this.TheoreticalMz = theoreticalMz;
            this.Identification = identification;
            this.RelativeMzError = (experimentalMz - theoreticalMz) / theoreticalMz;
        }
    }
}