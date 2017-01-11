namespace InternalLogic
{
    public class DataPoint
    {
        public double mz;
        public double rt;
        public int msnOrder;
        public double intensity;
        public int SelectedIonGuessChargeStateGuess;
        public double IsolationMZ;
        public double TotalIonCurrent;
        public double InjectionTime;
        public double relativeMZ;

        public DataPoint(double mz, double rt, int msnOrder, double intensity, double TotalIonCurrent, double InjectionTime, int SelectedIonGuessChargeStateGuess = 0, double IsolationMZ = 0, double relativeMZ = 0)
        {
            this.mz = mz;
            this.rt = rt;
            this.msnOrder = msnOrder;
            this.intensity = intensity;
            this.SelectedIonGuessChargeStateGuess = SelectedIonGuessChargeStateGuess;
            this.IsolationMZ = IsolationMZ;
            this.TotalIonCurrent = TotalIonCurrent;
            this.InjectionTime = InjectionTime;
            this.relativeMZ = relativeMZ;
        }

        public override string ToString()
        {
            return "(" + mz + "," + rt + ")";
        }
    }
}