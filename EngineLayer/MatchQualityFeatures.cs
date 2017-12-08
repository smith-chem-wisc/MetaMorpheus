namespace EngineLayer
{
    public class MatchQualityFeatures
    {
        #region Public Fields

        public readonly int matchingProductsHere;
        public readonly double intensityFracMatch;

        #endregion Public Fields

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