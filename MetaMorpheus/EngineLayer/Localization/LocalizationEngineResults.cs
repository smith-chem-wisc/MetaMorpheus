using System.Text;

namespace EngineLayer.Localization
{
    internal class LocalizationEngineResults : MetaMorpheusEngineResults
    {
        public LocalizationEngineResults(MetaMorpheusEngine s) : base(s)
        {
        }

        public LocalizationEngineResults(MetaMorpheusEngine s, int matchesWithTooManyPlacements) : base(s)
        {
            MatchesWithTooManyPlacements = matchesWithTooManyPlacements;
        }

        /// <summary>
        /// Matches whose modifications have more arrangements than the engine enumerates, and so have no placement evidence.
        /// </summary>
        public int MatchesWithTooManyPlacements { get; }

        public override string ToString()
        {
            var sb = new StringBuilder(base.ToString());
            if (MatchesWithTooManyPlacements > 0)
            {
                sb.AppendLine();
                sb.Append("Matches with too many modification arrangements to localize: " + MatchesWithTooManyPlacements);
            }
            return sb.ToString();
        }
    }
}
