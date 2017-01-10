using System;

namespace OldInternalLogic
{
    public class MSPeak
    {
        public double MZ { get; private set; }

        public double Intensity { get; private set; }

        public int Charge { get; private set; }

        public double mz { get; private set; }

        public MSPeak(double mz, double intensity, int charge, int polarity)
        {
            MZ = mz;
            Intensity = intensity;
            Charge = charge;
            CalculateMass(charge, polarity);
        }

        private void CalculateMass(int charge, int polarity)
        {
            if (charge == 0)
            {
                charge = polarity;
            }
            mz = MassFromMZ(MZ, charge);
        }

        public static double MassFromMZ(double mz, int charge)
        {
            return charge == 0 ? mz : mz * Math.Abs(charge) - charge * Constants.PROTON_MASS;
        }

        public static double MZFromMass(double mass, int charge)
        {
            if (charge == 0)
            {
                throw new ArgumentOutOfRangeException("Charge cannot be zero.");
            }
            return (mass + charge * Constants.PROTON_MASS) / Math.Abs(charge);
        }

        public static int AscendingMZComparison(MSPeak left, MSPeak right)
        {
            return left.MZ.CompareTo(right.MZ);
        }

        public static int DescendingIntensityComparison(MSPeak left, MSPeak right)
        {
            return -(left.Intensity.CompareTo(right.Intensity));
        }

        public static int AscendingMassComparison(MSPeak left, MSPeak right)
        {
            return left.mz.CompareTo(right.mz);
        }
    }
}