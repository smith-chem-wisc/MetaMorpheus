namespace EngineLayer.CrosslinkSearch
{
    class BestPeptideScoreNotch
    {
        public CompactPeptide BestPeptide { get; set; }
        public double BestScore { get; set; }
        public int BestNotch { get; set; }
        public BestPeptideScoreNotch(CompactPeptide bestPeptide, double bestScore, int bestNotch)
        {
            BestPeptide = bestPeptide;
            BestScore = bestScore;
            BestNotch = bestNotch;
        }
        public int[] topPosition { get; set; }
    }
}
