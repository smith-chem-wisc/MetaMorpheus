using Chemistry;
using MassSpectrometry;

namespace EngineLayer
{
    public class MatchedFragmentIon
    {
        public readonly TheoreticalFragmentIon TheoreticalFragmentIon;
        public readonly double Mz;
        public readonly double Intensity;
        public readonly double PpmMassError;

        /// <summary>
        /// Constructs a new MatchedFragmentIon given information about a theoretical and an experimental fragment mass spectral peak
        /// </summary>
        public MatchedFragmentIon(TheoreticalFragmentIon theoreticalFragmentIon, double experMz, double experIntensity)
        {
            this.TheoreticalFragmentIon = theoreticalFragmentIon;
            this.Mz = experMz;
            this.Intensity = experIntensity;
            this.PpmMassError = ((experMz.ToMass(theoreticalFragmentIon.Charge) - theoreticalFragmentIon.Mass) / theoreticalFragmentIon.Mass) * 1e6;
        }

        public int IntensityRank { get; set; }
        /// <summary>
        /// Summarizes a TheoreticalFragmentIon into a string for debug purposes
        /// TODO: Convert to a usable format for output
        /// </summary>
        public override string ToString()
        {
            return TheoreticalFragmentIon.ProductType.ToString().ToLowerInvariant() + TheoreticalFragmentIon.IonNumber + "+" + TheoreticalFragmentIon.Charge + "\t;" + TheoreticalFragmentIon.Mass;
        }
    }
}