using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Gptmd
{
    public class GptmdResults : MetaMorpheusEngineResults
    {
        private readonly int ModsAdded;

        public GptmdResults(MetaMorpheusEngine s, Dictionary<string, HashSet<Tuple<int, Modification>>> mods, int modsAdded) : base(s)
        {
            Mods = mods;
            ModsAdded = modsAdded;
        }

        public Dictionary<string, HashSet<Tuple<int, Modification>>> Mods { get; private set; }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("Modifications trying to add: " + ModsAdded);
            sb.AppendLine("Proteins trying to expand: " + Mods.Count);
            sb.AppendLine("Mods types and counts:");
            sb.Append(string.Join(Environment.NewLine, Mods.SelectMany(b => b.Value).GroupBy(b => b.Item2).OrderBy(b => -b.Count()).Select(b => "\t" + b.Key.IdWithMotif + "\t" + b.Count())));
            return sb.ToString();
        }
    }
}