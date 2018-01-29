using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class GptmdParameters
    {
        #region Public Constructors

        public GptmdParameters()
        {
            ListOfModsGptmd = GlobalEngineLevelSettings.AllModsKnown.Where(b =>
                b.modificationType.Equals("N-linked glycosylation") ||
                b.modificationType.Equals("Other glycosylation") ||
                b.modificationType.Equals("Mod") ||
                b.modificationType.Equals("PeptideTermMod") ||
                b.modificationType.Equals("Metal") ||
                b.modificationType.Equals("ProteinTermMod")).Select(b => (b.modificationType, b.id)).ToList();
        }

        #endregion Public Constructors

        #region Public Properties

        public List<(string, string)> ListOfModsGptmd { get; set; }

        #endregion Public Properties
    }
}