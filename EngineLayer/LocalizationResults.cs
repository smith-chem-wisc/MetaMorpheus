using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class LocalizationResults
    {

        #region Public Constructors

        public LocalizationResults(Dictionary<ProductType, double[]> matchedIonDictPositiveIsMatch, List<double> localizedScores)
        {
            this.MatchedIonMassesListPositiveIsMatch = new MatchedIonMassesListPositiveIsMatch(matchedIonDictPositiveIsMatch);
            this.LocalizedScores = localizedScores;
        }

        #endregion Public Constructors

        #region Public Properties

        public MatchedIonMassesListPositiveIsMatch MatchedIonMassesListPositiveIsMatch { get; }
        public List<double> LocalizedScores { get; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append(string.Join(";", MatchedIonMassesListPositiveIsMatch.Select(b => b.Value.Count(c => c > 0))) + '\t');

            sb.Append("[");
            foreach (var kvp in MatchedIonMassesListPositiveIsMatch)
                sb.Append("[" + string.Join(",", kvp.Value.Where(b => b > 0).Select(b => b.ToString("F5", CultureInfo.InvariantCulture))) + "];");
            sb.Append("]" + '\t');

            sb.Append("[" + string.Join(",", LocalizedScores.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))) + "]" + '\t');

            return sb.ToString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal static string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();
            sb.Append("Matched Ion Counts" + '\t');
            sb.Append("Matched Ion Masses" + '\t');
            sb.Append("Localized Scores");
            return sb.ToString();
        }

        #endregion Internal Methods

    }
}