using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class MatchedIonMassesListOnlyMatches : IEnumerable<KeyValuePair<ProductType, double[]>>, IEquatable<MatchedIonMassesListOnlyMatches>
    {
        #region Private Fields

        private readonly Dictionary<ProductType, double[]> matchedIonDictPositiveIsMatch;

        #endregion Private Fields

        #region Public Constructors

        public MatchedIonMassesListOnlyMatches(Dictionary<ProductType, double[]> matchedIonDictPositiveIsMatch)
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

        public override int GetHashCode()
        {
            return matchedIonDictPositiveIsMatch.SelectMany(b => b.Value).Count(b => b > 0);
        }

        public bool Equals(MatchedIonMassesListOnlyMatches other)
        {
            foreach (var hah in other)
            {
                if (!matchedIonDictPositiveIsMatch.TryGetValue(hah.Key, out double[] val))
                    return false;
                foreach (var ok in hah.Value.Where(b => b > 0))
                    if (!val.Contains(ok))
                        return false;
            }
            return true;
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