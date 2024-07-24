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
        public IBioPolymerWithSetMods Peptide { get; }
        public List<SpectralMatch> SpectralMatches { get; }

        /// <summary>
        /// This class groups all spectra associated with a given peptide together,
        /// to facilitate the calculation of PEP values.
        /// </summary>
        /// <param name="peptide"></param>
        /// <param name="spectralMatches"></param>
        public PeptideMatchGroup(IBioPolymerWithSetMods peptide, List<SpectralMatch> spectralMatches)
        {
            Peptide = peptide;
            SpectralMatches = spectralMatches;
        }

        public SpectralMatch BestMatch => SpectralMatches.MaxBy(match => match.Score);

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
