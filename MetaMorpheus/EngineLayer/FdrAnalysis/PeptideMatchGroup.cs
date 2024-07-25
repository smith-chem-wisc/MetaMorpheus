using Omics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class PeptideMatchGroup : IEnumerable<SpectralMatch>
    {
        public string PeptideFullSequence { get; }
        public List<SpectralMatch> SpectralMatches { get; }

        /// <summary>
        /// This class groups all spectra associated with a given peptide together,
        /// to facilitate the calculation of PEP values.
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="spectralMatches"></param>
        public PeptideMatchGroup(string fullPeptideSeq, List<SpectralMatch> spectralMatches)
        {
            PeptideFullSequence = fullPeptideSeq;
            SpectralMatches = spectralMatches;
        }

        public int TargetCount => SpectralMatches.Sum(p => p.BestMatchingBioPolymersWithSetMods.Count(t => !t.Peptide.Parent.IsDecoy));
        public int DecoyCount => SpectralMatches.Sum(p => p.BestMatchingBioPolymersWithSetMods.Count(t => t.Peptide.Parent.IsDecoy));

        public SpectralMatch BestMatch => SpectralMatches.MaxBy(match => match);

        public SpectralMatch BestMatchByPep => SpectralMatches.MinBy(match => match.FdrInfo.PEP);

        public IEnumerator<SpectralMatch> GetEnumerator()
        {
            return SpectralMatches.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }
}
