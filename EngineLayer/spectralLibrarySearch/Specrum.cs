using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;

namespace EngineLayer.spectralLibrarySearch
{
    public class Spectrum
    {
        public Spectrum()
        {
        }
        public Spectrum(String sequenceWithCharge, List<MatchedFragmentIon> peaks)
        {
            SequenceWithCharge = sequenceWithCharge; ;
            MatchedFragmentIons = peaks;
        }
        public Spectrum(String sequenceWithCharge,String sequence, double precursorMz, int charge_state, List<MatchedFragmentIon> peaks, Boolean isDecoy)
        {
            SequenceWithCharge = sequenceWithCharge;
            Sequence = sequence;
            PrecursorMz = precursorMz;
            MatchedFragmentIons = peaks;
            Charge_state = charge_state;
            IsDecoy = isDecoy;
        }

        public string SequenceWithCharge { get; set; }
        public string Sequence { get; set; }
        public double MW { get; set; }
        public double rententionTime { get; set; }
        public double PrecursorMz { get; set; }
        public int Charge_state { get; set; }
        public double totalIonCurrent { get; set; }
        public double MonoisotopicMass { get; set; }
        public List<MatchedFragmentIon> MatchedFragmentIons { get; set; }
        public Boolean IsDecoy { get; set; }


       
        public override string ToString()
        {
            StringBuilder spectrum = new StringBuilder();
            spectrum.Append("Name: " + Sequence + "/" + Charge_state + "\r\n");
            spectrum.Append("precursor: " + PrecursorMz);
            spectrum.Append("Matched peaks number : " + this.MatchedFragmentIons.Count + "\r\n");
            var intensitySum = this.MatchedFragmentIons.Select(b => b.Intensity).Sum();
            foreach (var eachPeak in this.MatchedFragmentIons)
            {
                var norIntensity = eachPeak.Intensity / intensitySum;
                spectrum.Append(eachPeak.Mz + "\t" + norIntensity + "\t" + "\"" + eachPeak.NeutralTheoreticalProduct.ProductType.ToString() + eachPeak.NeutralTheoreticalProduct.FragmentNumber.ToString() +"^"+ eachPeak.Charge+"/" + 0 + "ppm" + "\"" + "\r\n");
            }

            return spectrum.ToString();

        }
    }
}
