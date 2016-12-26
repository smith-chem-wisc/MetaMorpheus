namespace mzCal
{
    public class LabeledDataPoint
    {
        public double[] inputs;
        public double output;

        public LabeledDataPoint(double[] v1, double v2)
        {
            inputs = v1;
            output = v2;
        }
    }
}