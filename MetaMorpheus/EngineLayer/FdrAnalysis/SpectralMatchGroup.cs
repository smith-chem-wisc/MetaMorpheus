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
    public class SpectralMatchGroup : IEnumerable<SpectralMatch>
    {
        public string BaseSequence { get; }
        public List<SpectralMatch> SpectralMatches { get; }

        /// <summary>
        /// This class groups all spectral matches associated with a given sequence
        /// to facilitate the calculation of PEP values.
        /// </summary>
        /// <param name="baseSequence"> The sequence to be used for grouping</param>
        /// <param name="spectralMatches"> Every spectral match that matches the sequence</param>
        public SpectralMatchGroup(string baseSequence, List<SpectralMatch> spectralMatches)
        {
            BaseSequence = baseSequence;
            SpectralMatches = spectralMatches;
        }

        public static List<SpectralMatchGroup> GroupByBaseSequence(List<SpectralMatch> spectralMatches)
        {
            // This groups psms by base sequence, ensuring that PSMs with the same base sequence but different modifications are grouped together when training.

            // TODO: Determine if it's better to group PSMs by base sequence or by full sequence.
            return spectralMatches.GroupBy(p => p.BaseSequence)
                .Select(group => new SpectralMatchGroup(group.Key, group.ToList()))
                .OrderByDescending(matchGroup => matchGroup.Count())
                .ThenByDescending(matchGroup => matchGroup.BestMatch.Score)
                .ToList();
        }

        /// <summary>
        /// Returns the top-scoring PSM for each full sequence in the spectral match group.
        /// i.e., if the same base sequence has multiple modifications, this function will return the PSM with the highest score for each modification.
        /// </summary>
        public IEnumerable<SpectralMatch> GetBestMatches()
        {
            return SpectralMatches
                .GroupBy(p => p.FullSequence)
                .Select(g => g.MaxBy(p => p));
        }

        /// <summary>
        /// One match per peptide: each full sequence's best match by the default comparer, which reads the score and
        /// never PEP. It is the match <see cref="GetBestMatches"/> gives PEP to train on, so every peptide-level step
        /// after PEP stands each peptide on the match the model learned from, whatever PEPs its other matches got.
        /// Returned in descending order of those matches, the order the classic q-value pass needs.
        /// </summary>
        public static List<SpectralMatch> PeptideStandIns(IEnumerable<SpectralMatch> spectralMatches) =>
            PeptideStandIns(spectralMatches, p => p.FullSequence);

        /// <summary>
        /// As <see cref="PeptideStandIns(IEnumerable{SpectralMatch})"/>, with <paramref name="peptideOf"/> saying which
        /// matches are one peptide.
        /// </summary>
        public static List<SpectralMatch> PeptideStandIns<TKey>(IEnumerable<SpectralMatch> spectralMatches, Func<SpectralMatch, TKey> peptideOf) =>
            spectralMatches
                .OrderByDescending(p => p)
                .GroupBy(peptideOf)
                .Select(g => g.First()) // OrderByDescending is stable, so ties keep their order, as MaxBy's first maximum does
                .ToList();

        /// <summary>
        /// This function is called if there aren't enough peptides to train at the peptide level
        /// </summary>
        /// <param name="spectralMatches"></param>
        /// <returns></returns>
        public static List<SpectralMatchGroup> GroupByIndividualPsm(List<SpectralMatch> spectralMatches)
        {
            return spectralMatches.Select(psm => new SpectralMatchGroup(psm.FullSequence, new List<SpectralMatch> { psm }))
                .ToList();
        }

        public SpectralMatch BestMatch => SpectralMatches.MaxBy(match => match);

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