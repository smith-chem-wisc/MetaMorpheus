namespace MetaMorpheus
{
    public class MassTolerance
    {
        public double Value { get; private set; }

        public MassToleranceUnits Units { get; private set; }

        public MassTolerance(double value, MassToleranceUnits units)
        {
            Value = value;
            Units = units;
        }

        public static double CalculateMassError(double experimental, double theoretical, MassToleranceUnits massErrorUnits)
        {
            if (massErrorUnits == MassToleranceUnits.Da)
            {
                return experimental - theoretical;
            }
            else if (massErrorUnits == MassToleranceUnits.ppm)
            {
                return (experimental - theoretical) / theoretical * 1e6;
            }
            else
            {
                return double.NaN;
            }
        }

        public static double operator +(double left, MassTolerance right)
        {
            if (right.Units == MassToleranceUnits.Da)
            {
                return left + right.Value;
            }
            else if (right.Units == MassToleranceUnits.ppm)
            {
                return left + left * right.Value / 1e6;
            }
            else
            {
                return double.NaN;
            }
        }

        public static double operator -(double left, MassTolerance right)
        {
            if (right.Units == MassToleranceUnits.Da)
            {
                return left - right.Value;
            }
            else if (right.Units == MassToleranceUnits.ppm)
            {
                return left - left * right.Value / 1e6;
            }
            else
            {
                return double.NaN;
            }
        }
    }
}