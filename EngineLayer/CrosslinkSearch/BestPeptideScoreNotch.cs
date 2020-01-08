using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.CrosslinkSearch
{
    internal class BestPeptideScoreNotch
    {
        public BestPeptideScoreNotch(PeptideWithSetModifications bestPeptide, double bestScore, int bestNotch)
        {
            BestPeptide = bestPeptide;
            BestScore = bestScore;
            BestNotch = bestNotch;
        }

        public PeptideWithSetModifications BestPeptide { get; set; }
        public double BestScore { get; set; }
        public int BestNotch { get; set; }
        public int[] TopPosition { get; set; }
    }
}