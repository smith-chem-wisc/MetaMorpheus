namespace EngineLayer.Calibration
{
    public class LabeledMs2DataPoint : IHasInputsAndOutputs
    {

        #region Public Fields

        public readonly double mz;
        public readonly double rt;
        public readonly double intensity;
        public readonly double totalIonCurrent;
        public readonly double injectionTime;
        public readonly double isolationMz;
        public readonly Psm identification;

        #endregion Public Fields

        #region Public Constructors

        public LabeledMs2DataPoint(double mz, double rt, double intensity, double totalIonCurrent, double? injectionTime, double isolationMz, double label, Psm identification)
        {
            this.mz = mz;
            this.rt = rt;
            this.intensity = intensity;
            this.totalIonCurrent = totalIonCurrent;
            this.injectionTime = injectionTime ?? double.NaN;
            this.isolationMz = isolationMz;
            this.Label = label;
            this.identification = identification;
            Inputs = new double[] { mz, rt, intensity, totalIonCurrent, this.injectionTime, isolationMz };
        }

        #endregion Public Constructors

        #region Public Properties

        public static string TabSeparatedHeaderForMs1 { get { return "mz\trt\tintensity\tTIC\tInjectionTime\tIsolationmz\tLabel"; } }
        public double Label { get; private set; }
        public double[] Inputs { get; private set; }

        #endregion Public Properties

    }
}