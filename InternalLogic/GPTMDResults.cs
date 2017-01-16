using System;
using System.Collections.Generic;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class GPTMDResults : MyResults
    {
        #region Private Fields

        private readonly int modsAdded;

        #endregion Private Fields

        #region Public Constructors

        public GPTMDResults(MyEngine s, Dictionary<string, HashSet<Tuple<int, string>>> mods, int modsAdded) : base(s)
        {
            this.mods = mods;
            this.modsAdded = modsAdded;
        }

        #endregion Public Constructors

        #region Public Properties

        public Dictionary<string, HashSet<Tuple<int, string>>> mods { get; private set; }

        #endregion Public Properties

        #region Protected Methods

        protected override string GetStringForOutput()
        {
            var sb = new StringBuilder();
            sb.AppendLine("\t\tModifications added = " + modsAdded);
            sb.Append("\t\tProteins expanded = " + mods.Count);
            return sb.ToString();
        }

        #endregion Protected Methods
    }
}