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

        public static List<PeptideMatchGroup> GroupByFullSequence(List<SpectralMatch> spectralMatches)
        {
            // This groups psms by full sequences. All ambiguous PSMs are grouped together (full sequence == null), this prevents accidental cross contamination during training.
            return spectralMatches.GroupBy(p => p.FullSequence)
                .Select(group => new PeptideMatchGroup(group.Key, group.ToList()))
                .OrderByDescending(matchGroup => matchGroup.Count())
                .ThenByDescending(matchGroup => matchGroup.BestMatch.Score)
                .ToList();
        }

        /// <summary>
        /// This function is called if there aren't enough peptides to train at the peptide level
        /// </summary>
        /// <param name="spectralMatches"></param>
        /// <returns></returns>
        public static List<PeptideMatchGroup> GroupByIndividualPsm(List<SpectralMatch> spectralMatches)
        {
            return spectralMatches.Select(psm => new PeptideMatchGroup(psm.FullSequence, new List<SpectralMatch> { psm }))
                .ToList();
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