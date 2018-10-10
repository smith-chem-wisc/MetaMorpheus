using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.ModificationAnalysis
{
    public class ModificationAnalysisResults : MetaMorpheusEngineResults
    {
        public ModificationAnalysisResults(ModificationAnalysisEngine modificationAnalysisEngine) : base(modificationAnalysisEngine)
        {
        }

        /// <summary>
        /// String is the mod ID, integer is the count of that mod observed
        /// </summary>
        public Dictionary<string, int> CountOfAmbiguousButLocalizedModsSeen { get; internal set; }
        public Dictionary<string, int> CountOfModsSeenAndLocalized { get; internal set; }
        public Dictionary<string, int> CountOfEachModSeenOnProteins { get; internal set; }
        public Dictionary<string, int> CountOfUnlocalizedMods { get; internal set; }
        public Dictionary<ChemicalFormula, int> CountOfUnlocalizedFormulas { get; internal set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("Localized mods seen below q-value 0.01:");
            sb.AppendLine(string.Join(Environment.NewLine, CountOfModsSeenAndLocalized.OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            sb.AppendLine("(Approx) Additional localized but protein ambiguous mods seen below q-value 0.01:");
            sb.AppendLine(string.Join(Environment.NewLine, CountOfAmbiguousButLocalizedModsSeen.OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            sb.AppendLine("(Approx) Additional unlocalized mods seen below q-value 0.01:");
            sb.AppendLine(string.Join(Environment.NewLine, CountOfUnlocalizedMods.OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            sb.AppendLine("(Approx) Additional unlocalized modification formulas seen below q-value 0.01:");
            sb.AppendLine(string.Join(Environment.NewLine, CountOfUnlocalizedFormulas.OrderBy(b => -b.Value).Select(b => "\t" + b.Key.Formula + "\t" + b.Value)));
            sb.AppendLine();
            sb.AppendLine("All mods in database limited to peptides observed in the results:");
            sb.AppendLine(string.Join(Environment.NewLine, CountOfEachModSeenOnProteins.OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));

            return sb.ToString();
        }
    }
}