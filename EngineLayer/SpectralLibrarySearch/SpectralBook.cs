using Proteomics.Fragmentation;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.SpectralLibrarySearch
{
    //A spectral library is composed of books, where each book contains the MS2 spectrum and the ID
    public class SpectralBook
    {
        public SpectralBook(string sequence, double precursorMz, int charge_state, List<MatchedFragmentIon> peaks, bool isDecoy=false)
        {
            Sequence = sequence;
            PrecursorMz = precursorMz;
            MatchedFragmentIons = peaks;
            ChargeState = charge_state;
            IsDecoy = isDecoy;
            SequenceWithCharge = sequence + charge_state.ToString();
        }

        public string SequenceWithCharge { get; set; }
        public string Sequence { get; set; }
        public double RetentionTime { get; set; }
        public double PrecursorMz { get; set; }
        public int ChargeState { get; set; }
        public double TotalIonCurrent { get; set; }
        public double MonoisotopicMass { get; set; }
        public List<MatchedFragmentIon> MatchedFragmentIons { get; set; }
        public bool IsDecoy { get; set; }



        public override string ToString()
        {
            StringBuilder spectrum = new StringBuilder();
            spectrum.Append("Name: " + Sequence + "/" + ChargeState + "\r\n");
            spectrum.Append("precursor: " + PrecursorMz);
            spectrum.Append("Matched peaks number : " + this.MatchedFragmentIons.Count + "\r\n");
            var intensitySum = MatchedFragmentIons.Select(b => b.Intensity).Sum();
            foreach (var eachPeak in this.MatchedFragmentIons)
            {
                var norIntensity = eachPeak.Intensity / intensitySum;
                spectrum.Append(eachPeak.Mz + "\t" + norIntensity + "\t" + "\"" + eachPeak.NeutralTheoreticalProduct.ProductType.ToString() + eachPeak.NeutralTheoreticalProduct.FragmentNumber.ToString() + "^" + eachPeak.Charge + "/" + 0 + "ppm" + "\"" + "\r\n");
            }

            return spectrum.ToString();
        }
    }
}