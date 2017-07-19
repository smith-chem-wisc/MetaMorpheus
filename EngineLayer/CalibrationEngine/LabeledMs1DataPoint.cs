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
        public readonly SingleScanManyPeptidesMatch identification;

        #endregion Public Fields

        #region Public Constructors

        public LabeledMs1DataPoint(double mz, double rt, double intensity, double totalIonCurrent, double? injectionTime, double label, SingleScanManyPeptidesMatch identification)
        {
            this.mz = mz;
            this.rt = rt;
            this.intensity = intensity;
            this.totalIonCurrent = totalIonCurrent;
            this.injectionTime = injectionTime ?? double.NaN;
            this.Label = label;
            this.identification = identification;
            Inputs = new double[] { mz, rt, intensity, totalIonCurrent, this.injectionTime };
        }

        #endregion Public Constructors

        #region Public Properties

        public static string TabSeparatedHeaderForMs1 { get { return "mz\trt\tintensity\tTIC\tInjectionTime\tLabel"; } }
        public double Label { get; private set; }
        public double[] Inputs { get; private set; }

        #endregion Public Properties

    }
}