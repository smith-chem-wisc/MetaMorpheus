namespace EngineLayer.Calibration
{
    public class LabeledMs2DataPoint : IHasInputsAndOutputs
    {
        #region Public Fields

        public readonly double mz;
        public readonly double rt;
        
        public readonly double LOGinjectionTime;

        public readonly double LOGtotalIonCurrent;
        public readonly Psm identification;

        #endregion Public Fields

        #region Public Constructors

        public LabeledMs2DataPoint(double mz, double rt, double LOGtotalIonCurrent, double LOGinjectionTime, double label, Psm identification)
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

        public static string TabSeparatedHeaderForMs1 { get { return "mz\trt\tLOGTIC\tLOGInjectionTime\tLabel"; } }

        public double Label { get; private set; }
        public double[] Inputs { get; private set; }

        #endregion Public Properties
    }
}