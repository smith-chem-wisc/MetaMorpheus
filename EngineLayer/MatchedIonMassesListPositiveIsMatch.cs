using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class MatchedIonMassesListPositiveIsMatch : IEnumerable<KeyValuePair<ProductType, double[]>>
    {

        #region Private Fields

        private Dictionary<ProductType, double[]> matchedIonDictPositiveIsMatch;

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
                double[] val;
                if (!matchedIonDictPositiveIsMatch.TryGetValue(hah.Key, out val))
                    return false;
                for (int i = 0; i < val.Length; i++)
                {
                    if (val[i] != hah.Value[i])
                        return false;
                }
            }
            return true;
        }

        public override int GetHashCode()
        {
            return matchedIonDictPositiveIsMatch.SelectMany(b => b.Value).Count();
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