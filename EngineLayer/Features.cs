namespace EngineLayer
{
    public class MatchQualityFeatures
    {
        #region Private Fields

        private readonly int matchingProductsHere;
        private readonly double intensityFracMatch;

        #endregion Private Fields

        #region Public Constructors

        public MatchQualityFeatures(int matchingProductsHere, double intensityFracMatch)
        {
            this.matchingProductsHere = matchingProductsHere;
            this.intensityFracMatch = intensityFracMatch;
        }

        public MatchQualityFeatures()
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return matchingProductsHere + "\t" + intensityFracMatch;
        }

        public double[] ToDoubleArray()
        {
            return new[] { matchingProductsHere, intensityFracMatch };
        }

        #endregion Public Methods
    }
}