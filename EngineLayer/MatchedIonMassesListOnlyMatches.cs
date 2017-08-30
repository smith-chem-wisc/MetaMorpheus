using System.Collections;
using System.Collections.Generic;

namespace EngineLayer
{
    public class MatchedIonMassesListOnlyMasses : IEnumerable<KeyValuePair<ProductType, double[]>>

    {
        #region Private Fields

        private readonly Dictionary<ProductType, double[]> matchedIonDictPositiveIsMatch;

        #endregion Private Fields

        #region Public Constructors

        public MatchedIonMassesListOnlyMasses(Dictionary<ProductType, double[]> matchedIonDictPositiveIsMatch)
        {
            this.matchedIonDictPositiveIsMatch = matchedIonDictPositiveIsMatch;
        }

        #endregion Public Constructors

        #region Public Indexers

        public double[] this[ProductType y]
        {
            get
            {
                return matchedIonDictPositiveIsMatch[y];
            }
        }

        #endregion Public Indexers

        #region Public Methods

        public IEnumerator<KeyValuePair<ProductType, double[]>> GetEnumerator()
        {
            return matchedIonDictPositiveIsMatch.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return matchedIonDictPositiveIsMatch.GetEnumerator();
        }

        #endregion Public Methods

        #region Internal Methods

        internal bool ContainsKey(ProductType y)
        {
            return matchedIonDictPositiveIsMatch.ContainsKey(y);
        }

        #endregion Internal Methods
    }
}