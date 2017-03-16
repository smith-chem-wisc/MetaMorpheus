using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.Gptmd
{
    public class GptmdResults : MyResults
    {

        #region Private Fields

        private readonly int modsAdded;

        #endregion Private Fields

        #region Public Constructors

        public GptmdResults(MyEngine s, Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>> mods, int modsAdded) : base(s)
        {
            this.Mods = mods;
            this.modsAdded = modsAdded;
        }

        #endregion Public Constructors

        #region Public Properties

        public Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>> Mods { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("Modifications added = " + modsAdded);
            sb.AppendLine("Proteins expanded = " + Mods.Count);
            sb.AppendLine("Mods added types and counts:");
            sb.Append(string.Join(Environment.NewLine, Mods.SelectMany(b => b.Value).GroupBy(b => b.Item2).OrderBy(b => -b.Count()).Select(b => b.Key.id + " " + b.Count())));
            return sb.ToString();
        }

        #endregion Public Methods

    }
}