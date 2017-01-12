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

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("Modifications added = " + modsAdded);
            sb.AppendLine("Proteins expanded = " + mods.Count);
            return sb.ToString();
        }
    }
}