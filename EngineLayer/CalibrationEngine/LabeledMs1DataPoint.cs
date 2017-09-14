namespace EngineLayer.Calibration
{
    public class LabeledMs1DataPoint : IHasInputsAndOutputs
    {
        #region Public Fields

        public readonly double mz;
        public readonly double rt;
        
        public readonly double LOGtotalIonCurrent;

        public readonly double LOGinjectionTime;
        public readonly Psm identification;

        #endregion Public Fields

        #region Public Constructors

        public LabeledMs1DataPoint(double mz, double rt, double LOGtotalIonCurrent, double LOGinjectionTime, double label, Psm identification)
        {
            this.mz = mz;
            this.rt = rt;
            this.LOGtotalIonCurrent = LOGtotalIonCurrent;
            this.LOGinjectionTime = LOGinjectionTime;
            this.Label = label;
            this.identification = identification;
            Inputs = new double[] { this.mz, this.rt, this.LOGtotalIonCurrent, this.LOGinjectionTime };
        }

        #endregion Public Constructors

        #region Public Properties

        //public static string TabSeparatedHeaderForMs1 { get { return "mz\trt\tintensity\tTIC\tInjectionTime\tLabel"; } }
        public static string TabSeparatedHeaderForMs1 { get { return "mz\trt\tLOGTIC\tLOGinjectionTime\tLabel"; } }

        public double Label { get; private set; }
        public double[] Inputs { get; private set; }

        #endregion Public Properties
    }
}