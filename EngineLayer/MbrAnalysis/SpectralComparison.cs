using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry.MzSpectra;

namespace EngineLayer.MbrAnalysis
{
    public class SpectralComparison
    {
        public string SpectraSource;
        public string LibrarySpectraSource;
        public string FullSequence;
        public Dictionary<string, double> ScoreDictionary;

        public SpectralComparison(string spectraSource, string librarySpectraSource, string fullSequence, SpectralSimilarity similarity)
        {
            SpectraSource = spectraSource;
            LibrarySpectraSource = librarySpectraSource;
            FullSequence = fullSequence;
            ScoreDictionary = MakeScoreDictionary(similarity);
        }

        private Dictionary<string, double> MakeScoreDictionary(SpectralSimilarity similarity)
        {
            ScoreDictionary = new Dictionary<string, double>();
            ScoreDictionary.Add("")
        }
    }
}
