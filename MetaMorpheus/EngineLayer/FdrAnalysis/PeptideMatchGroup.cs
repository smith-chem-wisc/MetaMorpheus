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
        /// <param name="fullPeptideSeq"> The full sequence to be used for grouping</param>
        /// <param name="spectralMatches"> Every spectral match that matches the full sequence</param>
        public PeptideMatchGroup(string fullPeptideSeq, List<SpectralMatch> spectralMatches)
        {
            PeptideFullSequence = fullPeptideSeq;
            SpectralMatches = spectralMatches;
        }

        /// <summary>
        /// Returns the number of full sequences that match to at least one target protein.
        /// </summary>
        public int TargetCount => SpectralMatches.Sum(p => p.BestMatchingBioPolymersWithSetMods
            .Select(t => t.Peptide)
            .GroupBy(peptide => peptide.FullSequence)
            .Count(group => group.Any(p => !p.Parent.IsDecoy)));

        /// <summary>
        /// Returns the number of full sequences that match to at least one decoy protein.
        /// </summary>
        public int DecoyCount => SpectralMatches.Sum(p => p.BestMatchingBioPolymersWithSetMods
            .Select(t => t.Peptide)
            .GroupBy(peptide => peptide.FullSequence)
            .Count(group => group.Any(p => p.Parent.IsDecoy)));

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