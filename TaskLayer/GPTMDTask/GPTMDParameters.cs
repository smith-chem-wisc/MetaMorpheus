using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class GptmdParameters
    {
        public GptmdParameters()
        {
            ListOfModsGptmd = GlobalVariables.AllModsKnown.Where(b =>
                b.modificationType.Equals("N-linked glycosylation") ||
                b.modificationType.Equals("Other glycosylation") ||
                b.modificationType.Equals("Mod") ||
                //b.modificationType.Equals("Artifact") ||
                //b.modificationType.Equals("Biological") ||
                b.modificationType.Equals("PeptideTermMod") ||
                b.modificationType.Equals("Metal") ||
                b.modificationType.Equals("ProteinTermMod")).Select(b => (b.modificationType, b.id)).ToList();
        }

        public List<(string, string)> ListOfModsGptmd { get; set; }
    }
}