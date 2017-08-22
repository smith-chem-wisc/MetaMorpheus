using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public class ModificationAnalysisResults : MetaMorpheusEngineResults
    {
        #region Public Constructors

        public ModificationAnalysisResults(ModificationAnalysisEngine modificationAnalysisEngine) : base(modificationAnalysisEngine)
        {
        }

        #endregion Public Constructors

        #region Public Properties

        public Dictionary<string, int>[] AmbiguousButLocalizedModsSeen { get; internal set; }
        public Dictionary<string, int>[] ModsSeenAndLocalized { get; internal set; }
        public Dictionary<string, int>[] AllModsOnProteins { get; internal set; }
        public Dictionary<string, int>[] UnlocalizedMods { get; internal set; }
        public Dictionary<ChemicalFormula, int>[] UnlocalizedFormulas { get; internal set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            for (int i = 0; i < AllModsOnProteins.Length; i++)
            {
                sb.AppendLine("Search mode index:" + i + ". All mods in database limited to peptides observed in the results:");
                sb.AppendLine(string.Join(Environment.NewLine, AllModsOnProteins[i].OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
                sb.AppendLine("Search mode index:" + i + ". Localized mods seen below q-value 0.01:");
                sb.AppendLine(string.Join(Environment.NewLine, ModsSeenAndLocalized[i].OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
                sb.AppendLine("Search mode index:" + i + ". (Approx) Additional localized but protein ambiguous mods seen below q-value 0.01:");
                sb.AppendLine(string.Join(Environment.NewLine, AmbiguousButLocalizedModsSeen[i].OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
                sb.AppendLine("Search mode index:" + i + ". (Approx) Additional unlocalized mods seen below q-value 0.01:");
                sb.AppendLine(string.Join(Environment.NewLine, UnlocalizedMods[i].OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
                sb.AppendLine("Search mode index:" + i + ". (Approx) Additional unlocalized modification formulas seen below q-value 0.01:");
                sb.AppendLine(string.Join(Environment.NewLine, UnlocalizedFormulas[i].OrderBy(b => -b.Value).Select(b => "\t" + b.Key.Formula + "\t" + b.Value)));
            }
            return sb.ToString();
        }

        #endregion Public Methods
    }
}