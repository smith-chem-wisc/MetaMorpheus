using Proteomics.Fragmentation;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    /// <summary>
    /// This class describes one element of a spectral library.
    /// </summary>
    public class LibrarySpectrum
    {
        public string SequenceWithCharge { get; set; }
        public string Sequence { get; set; }
        public double RetentionTime { get; set; }
        public double PrecursorMz { get; set; }
        public int ChargeState { get; set; }
        public List<MatchedFragmentIon> MatchedFragmentIons { get; set; }
        public bool IsDecoy { get; set; }

        public LibrarySpectrum(string sequence, double precursorMz, int chargeState, List<MatchedFragmentIon> peaks, double rt, bool isDecoy = false)
        {
            Sequence = sequence;
            PrecursorMz = precursorMz;
            MatchedFragmentIons = peaks;
            ChargeState = chargeState;
            IsDecoy = isDecoy;
            RetentionTime = rt;
            SequenceWithCharge = sequence + chargeState.ToString();
        }

        public override string ToString()
        {
            StringBuilder spectrum = new StringBuilder();
            spectrum.Append("Name: " + Sequence + "/" + ChargeState);
            spectrum.Append("\nMW: " + PrecursorMz);
            spectrum.Append("\nComment: ");
            spectrum.Append("Parent=" + PrecursorMz);
            spectrum.Append(" Retention Time=" + RetentionTime);
            spectrum.Append("\nMatched peaks number : " + MatchedFragmentIons.Count + "\n");

            double maxIntensity = MatchedFragmentIons.Select(b => b.Intensity).Max();

            foreach (MatchedFragmentIon matchedIon in MatchedFragmentIons)
            {
                double intensityFraction = matchedIon.Intensity / maxIntensity;

                spectrum.Append(matchedIon.Mz + "\t" + intensityFraction + "\t" + "\"" +
                    matchedIon.NeutralTheoreticalProduct.ProductType.ToString() +
                    matchedIon.NeutralTheoreticalProduct.FragmentNumber.ToString() + "^" +
                    matchedIon.Charge + "/" + 0 + "ppm" + "\"" + "\n");
            }

            return spectrum.ToString();
        }
    }
}