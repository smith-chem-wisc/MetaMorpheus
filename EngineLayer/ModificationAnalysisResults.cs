using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    internal class ModificationAnalysisResults : MetaMorpheusEngineResults
    {

        #region Public Constructors

        public ModificationAnalysisResults(ModificationAnalysisEngine modificationAnalysisEngine) : base(modificationAnalysisEngine)
        {
        }

        #endregion Public Constructors

        #region Public Properties

        public Dictionary<string, int>[] AllModsSeen { get; internal set; }
        public Dictionary<string, int>[] AllModsOnPeptides { get; internal set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            for (int i = 0; i < AllModsOnPeptides.Length; i++)
            {
                sb.AppendLine("Search mode " + i + " Mods seen:");
                sb.AppendLine(string.Join(Environment.NewLine, AllModsSeen[i].OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
                sb.AppendLine("Search mode " + i + " Mods in database:");
                sb.AppendLine(string.Join(Environment.NewLine, AllModsOnPeptides[i].OrderBy(b => -b.Value).Select(b => "\t" + b.Key + "\t" + b.Value)));
            }
            return sb.ToString();
        }

        #endregion Public Methods

    }
}