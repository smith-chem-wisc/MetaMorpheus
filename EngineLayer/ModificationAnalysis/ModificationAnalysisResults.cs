using Chemistry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.ModificationAnalysis
{
    public class ModificationAnalysisResults : MetaMorpheusEngineResults
    {
        public ModificationAnalysisResults(ModificationAnalysisEngine modificationAnalysisEngine)
            : base(modificationAnalysisEngine)
        {
        }

        public Dictionary<string, int> AmbiguousButLocalizedModsSeen { get; internal set; }
        public Dictionary<string, int> ModsSeenAndLocalized { get; internal set; }
        public Dictionary<string, int> AllModsOnProteins { get; internal set; }
        public Dictionary<string, int> UnlocalizedMods { get; internal set; }
        public Dictionary<ChemicalFormula, int> UnlocalizedFormulas { get; internal set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("Localized mods seen below q-value 0.01:");
            sb.AppendLine(string.Join(Environment.NewLine, ModsSeenAndLocalized.OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            sb.AppendLine("(Approx) Additional localized but protein ambiguous mods seen below q-value 0.01:");
            sb.AppendLine(string.Join(Environment.NewLine, AmbiguousButLocalizedModsSeen.OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            sb.AppendLine("(Approx) Additional unlocalized mods seen below q-value 0.01:");
            sb.AppendLine(string.Join(Environment.NewLine, UnlocalizedMods.OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            sb.AppendLine("(Approx) Additional unlocalized modification formulas seen below q-value 0.01:");
            sb.AppendLine(string.Join(Environment.NewLine, UnlocalizedFormulas.OrderBy(b => -b.Value).Select(b => "\t" + b.Key.Formula + "\t" + b.Value)));
            sb.AppendLine();
            sb.AppendLine("All mods in database limited to peptides observed in the results:");
            sb.AppendLine(string.Join(Environment.NewLine, AllModsOnProteins.OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));

            return sb.ToString();
        }
    }
}