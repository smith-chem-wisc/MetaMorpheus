namespace EngineLayer.Calibration
{
    public class LabeledMs1DataPoint : IHasInputsAndOutputs
    {

        #region Public Fields

        public readonly double mz;
        public readonly double rt;
        public readonly double intensity;
        public readonly double totalIonCurrent;
        public readonly double injectionTime;

        #endregion Public Fields

        #region Public Constructors

        public LabeledMs1DataPoint(double mz, double rt, double intensity, double totalIonCurrent, double? injectionTime, double label)
        {
            this.mz = mz;
            this.rt = rt;
            this.intensity = intensity;
            this.totalIonCurrent = totalIonCurrent;
            this.injectionTime = injectionTime.HasValue ? injectionTime.Value : double.NaN;
            this.label = label;
            inputs = new double[] { mz, rt, intensity, totalIonCurrent, this.injectionTime };
        }

        #endregion Public Constructors

        #region Public Properties

        public static string TabSeparatedHeaderForMs1 { get { return "mz\trt\tintensity\tTIC\tInjectionTime\tLabel"; } }
        public double label { get; private set; }
        public double[] inputs { get; private set; }

        #endregion Public Properties

    }
}