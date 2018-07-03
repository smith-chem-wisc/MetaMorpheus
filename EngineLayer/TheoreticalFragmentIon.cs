using Chemistry;

namespace EngineLayer
{
    public class TheoreticalFragmentIon
    {
        public readonly double Mass;
        public readonly int Charge;
        public readonly double Mz;
        public readonly double Intensity;
        public readonly ProductType ProductType;
        public readonly int IonNumber;

        /// <summary>
        /// Constructs a new TheoreticalFragmentIon given information about its theoretical properties
        /// </summary>
        public TheoreticalFragmentIon(double mass, double theorIntensity, int charge, ProductType productType, int ionNumber)
        {
            this.Mass = mass;
            this.Charge = charge;
            this.Intensity = theorIntensity;
            this.IonNumber = ionNumber;
            this.ProductType = productType;
            this.Mz = mass.ToMz(charge);
        }

        /// <summary>
        /// Summarizes a TheoreticalFragmentIon into a string for debug purposes
        /// </summary>
        public override string ToString()
        {
            return ProductType + "" + IonNumber + ";" + Mass;
        }
    }
}
