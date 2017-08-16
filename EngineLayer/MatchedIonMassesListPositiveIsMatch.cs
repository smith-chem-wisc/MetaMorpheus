using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class MatchedIonMassesListPositiveIsMatch : IEnumerable<KeyValuePair<ProductType, double[]>>
    {
        #region Private Fields

        private readonly Dictionary<ProductType, double[]> matchedIonDictPositiveIsMatch;

        #endregion Private Fields

        #region Public Constructors

        public MatchedIonMassesListPositiveIsMatch(Dictionary<ProductType, double[]> matchedIonDictPositiveIsMatch)
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

        public override bool Equals(object obj)
        {
            var kk = obj as MatchedIonMassesListPositiveIsMatch;
            if (kk == null)
                return false;
            foreach (var hah in kk)
            {
                if (!matchedIonDictPositiveIsMatch.TryGetValue(hah.Key, out double[] val))
                    return false;
                foreach (var ok in hah.Value.Where(b => b > 0))
                    if (!val.Contains(ok))
                        return false;
            }
            return true;
        }

        public override int GetHashCode()
        {
            return matchedIonDictPositiveIsMatch.SelectMany(b => b.Value).Count(b => b > 0);
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