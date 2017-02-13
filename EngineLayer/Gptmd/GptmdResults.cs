using Proteomics;
using System;
using System.Collections.Generic;
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

        #region Protected Properties

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("Modifications added = " + modsAdded);
            sb.Append("Proteins expanded = " + Mods.Count);
            return sb.ToString();
        }

        #endregion Protected Properties

    }
}