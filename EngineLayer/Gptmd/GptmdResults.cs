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

        public GptmdResults(MyEngine s, Dictionary<string, HashSet<Tuple<int, BaseModification>>> mods, int modsAdded) : base(s)
        {
            this.Mods = mods;
            this.modsAdded = modsAdded;
        }

        #endregion Public Constructors

        #region Public Properties

        public Dictionary<string, HashSet<Tuple<int, BaseModification>>> Mods { get; private set; }

        #endregion Public Properties

        #region Protected Properties

        protected override string StringForOutput
        {
            get
            {
                var sb = new StringBuilder();
                sb.AppendLine("\t\tModifications added = " + modsAdded);
                sb.Append("\t\tProteins expanded = " + Mods.Count);
                return sb.ToString();
            }
        }

        #endregion Protected Properties

    }
}