namespace EngineLayer
{
    public class MatchQualityFeatures
    {
        #region Public Fields

        public readonly double[] arr;

        #endregion Public Fields

        #region Public Constructors

        public MatchQualityFeatures(double[] arr)
        {
            this.arr = arr;
        }

        public MatchQualityFeatures()
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return string.Join("\t", arr);
        }

        #endregion Public Methods
    }
}