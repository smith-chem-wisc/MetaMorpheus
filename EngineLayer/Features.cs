namespace EngineLayer
{
    public class Features
    {
        #region Private Fields

        private (double, double) score;

        #endregion Private Fields

        #region Public Constructors

        public Features((double, double) score)
        {
            this.score = score;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return score.Item1 + "\t" + score.Item2;
        }

        public double[] ToDoubleArray()
        {
            return new []{ score.Item1, score.Item2 };
        }

        #endregion Public Methods
    }
}