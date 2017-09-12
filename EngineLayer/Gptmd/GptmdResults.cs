using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Gptmd
{
    public class GptmdResults : MetaMorpheusEngineResults
    {
        #region Private Fields

        private readonly int modsAdded;

        #endregion Private Fields

        #region Public Constructors

        public GptmdResults(MetaMorpheusEngine s, Dictionary<string, HashSet<Tuple<int, Modification>>> mods, int modsAdded) : base(s)
        {
            this.Mods = mods;
            this.modsAdded = modsAdded;
        }

        #endregion Public Constructors

        #region Public Properties

        public Dictionary<string, HashSet<Tuple<int, Modification>>> Mods { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("Modifications trying to add: " + modsAdded);
            sb.AppendLine("Proteins trying to expand: " + Mods.Count);
            sb.AppendLine("Mods types and counts:");
            sb.Append(string.Join(Environment.NewLine, Mods.SelectMany(b => b.Value).GroupBy(b => b.Item2).OrderBy(b => -b.Count()).Select(b => "\t" + b.Key.id + "\t" + b.Count())));
            return sb.ToString();
        }

        #endregion Public Methods
    }
}