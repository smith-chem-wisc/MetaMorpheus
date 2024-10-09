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

        public IEnumerable<SpectralMatch> GetMatches(string searchType)
        {
            switch (searchType)
            {
                case "top-down":
                    return SpectralMatches;
                default:
                    return SpectralMatches.GroupBy(p => p.FullSequence).Select(g => g.MaxBy(p => p));
            }
        }

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