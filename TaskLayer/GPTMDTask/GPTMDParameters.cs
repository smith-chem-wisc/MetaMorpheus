using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class GptmdParameters
    {
        public GptmdParameters()
        {
            ListOfModsGptmd = GlobalVariables.AllModsKnown.Where(b =>
                b.ModificationType.Equals("N-linked glycosylation") ||
                b.ModificationType.Equals("Other glycosylation") ||
                b.ModificationType.Equals("Mod") ||
                //b.ModificationType.Equals("Artifact") ||
                //b.ModificationType.Equals("Biological") ||
                b.ModificationType.Equals("PeptideTermMod") ||
                b.ModificationType.Equals("Metal") ||
                b.ModificationType.Equals("ProteinTermMod")).Select(b => (b.ModificationType, b.IdWithMotif)).ToList();
        }

        public List<(string, string)> ListOfModsGptmd { get; set; }
    }
}