using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.spectralLibrarySearch
{
    public class SpectralLibrarySearchResults
    {

        public SpectralLibrarySearchResults(Spectrum theExperimentalSpectrum, int indexOfExperimentalSpectrum, List<SpectralLibrarayMatch> spectralLibrarayMatchs)
        {
            TheExperimentalSpectrum = theExperimentalSpectrum;
            IndexOfExperimentalSpectrum = indexOfExperimentalSpectrum;
            SpectralLibrarayMatchs = spectralLibrarayMatchs;
        }
        public SpectralLibrarySearchResults(string peptide, double score, int index)
        {
            Peptide = peptide;
            Score = score;
            Index = index;
        }

        public Spectrum TheExperimentalSpectrum { get; set; }
        public int IndexOfExperimentalSpectrum { get; set; }
        public List<SpectralLibrarayMatch> SpectralLibrarayMatchs { get; set; }
        public string Peptide { get; set; }
        public double Score { get; set; }
        public int Index { get; set; }
    }
}
