using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.spectralLibrarySearch
{
    public class SpectralLibrarayMatch
    {

        //public Spectrum TheExperimentalSpectrum { get; set; }
        //public int IndexOfExperimentalSpectrum { get; set; }
        public Spectrum MatchedSpectrumFromLibrary { get; set; }
        public double MatchScore { get; set; }
        public List<MatchedFragmentIon> MatchedPeaks { get; set; }
        public int ScanIndex { get; set; }
        public CommonParameters CommonParameters { get; set; }
        public Ms2ScanWithSpecificMass Scan;
        public const double ToleranceForScoreDifferentiation = 1e-9;
        private List<(int Notch, Spectrum spectrum)> _BestLibraryMatchingPeptides;
        private double ScoreCutoff = 0.5;
        public double RunnerUpScore { get; set; }
        public Dictionary<Spectrum, List<MatchedFragmentIon>> SpectrumToMatchingFragments { get; private set; }
        public double Xcorr;

        public SpectralLibrarayMatch( Spectrum matchedSpectrumFromLibrary, int notch, double matchScore,int scanIndex, Ms2ScanWithSpecificMass scan, CommonParameters commonParameters, List<MatchedFragmentIon> matchedPeaks,double xcorr)
        {
            //TheExperimentalSpectrum = theExperimentalSpectrum;
            //IndexOfExperimentalSpectrum = indexOfExperimentalSpectrum;
            MatchedSpectrumFromLibrary = matchedSpectrumFromLibrary;
            MatchScore = matchScore;
            MatchedPeaks = matchedPeaks;
            ScanIndex = scanIndex;
            CommonParameters = commonParameters;
            Scan = scan;
            RunnerUpScore = ScoreCutoff;
            Xcorr = xcorr;
            AddOrReplace(matchedSpectrumFromLibrary, matchScore, notch, true, matchedPeaks, xcorr);

        }

        public void AddOrReplace(Spectrum matchedSpectrum, double newScore, int notch, bool reportAllAmbiguity, List<MatchedFragmentIon> matchedFragmentIons, double newXcorr)
        {
            if (newScore - MatchScore > ToleranceForScoreDifferentiation) //if new score beat the old score, overwrite it
            {
                _BestLibraryMatchingPeptides.Clear();
                _BestLibraryMatchingPeptides.Add((notch, matchedSpectrum));

                if (MatchScore - RunnerUpScore > ToleranceForScoreDifferentiation)
                {
                    RunnerUpScore = MatchScore;
                }

                MatchScore = newScore;
                Xcorr = newXcorr;

                SpectrumToMatchingFragments.Clear();
                SpectrumToMatchingFragments.Add(matchedSpectrum, matchedFragmentIons);
            }
            else if (newScore - MatchScore > -ToleranceForScoreDifferentiation && reportAllAmbiguity) //else if the same score and ambiguity is allowed
            {
                _BestLibraryMatchingPeptides.Add((notch, matchedSpectrum));

                if (!SpectrumToMatchingFragments.ContainsKey(matchedSpectrum))
                {
                    SpectrumToMatchingFragments.Add(matchedSpectrum, matchedFragmentIons);
                }
            }
            else if (newScore - RunnerUpScore > ToleranceForScoreDifferentiation)
            {
                RunnerUpScore = newScore;
            }
        }

       
    }
}
