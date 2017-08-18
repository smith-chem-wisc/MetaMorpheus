namespace EngineLayer.CrosslinkSearch
{
    internal class BestPeptideScoreNotch
    {
        #region Public Constructors

        public BestPeptideScoreNotch(CompactPeptide bestPeptide, double bestScore, int bestNotch)
        {
            BestPeptide = bestPeptide;
            BestScore = bestScore;
            BestNotch = bestNotch;
        }

        #endregion Public Constructors

        #region Public Properties

        public CompactPeptide BestPeptide { get; set; }
        public double BestScore { get; set; }
        public int BestNotch { get; set; }
        public int[] TopPosition { get; set; }

        #endregion Public Properties
    }
}