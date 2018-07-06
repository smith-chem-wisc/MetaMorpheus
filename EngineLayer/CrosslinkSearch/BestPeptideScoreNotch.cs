using Proteomics.ProteolyticDigestion;

namespace EngineLayer.CrosslinkSearch
{
    internal class BestPeptideScoreNotch
    {
        public BestPeptideScoreNotch(CompactPeptide bestPeptide, double bestScore, int bestNotch)
        {
            BestPeptide = bestPeptide;
            BestScore = bestScore;
            BestNotch = bestNotch;
        }

        public CompactPeptide BestPeptide { get; set; }
        public double BestScore { get; set; }
        public int BestNotch { get; set; }
        public int[] TopPosition { get; set; }
    }
}