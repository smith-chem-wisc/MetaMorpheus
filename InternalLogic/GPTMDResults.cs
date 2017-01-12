using System;
using System.Collections.Generic;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class GPTMDResults : MyResults
    {
        public Dictionary<string, HashSet<Tuple<int, string>>> mods { get; private set; }
        private int modsAdded;

        public GPTMDResults(MyEngine s, Dictionary<string, HashSet<Tuple<int, string>>> mods, int modsAdded) : base(s)
        {
            this.mods = mods;
            this.modsAdded = modsAdded;
        }

        protected override string GetStringForOutput()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("\t\tModifications added = " + modsAdded);
            sb.Append("\t\tProteins expanded = " + mods.Count);
            return sb.ToString();
        }
    }
}