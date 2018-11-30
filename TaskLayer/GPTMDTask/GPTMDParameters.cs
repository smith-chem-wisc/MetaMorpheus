using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class GptmdParameters
    {
        public GptmdParameters()
        {
            ListOfModsGptmd = GlobalVariables.AllModsKnown.Where(b =>
                b.ModificationType.Equals("Common Artifact") ||
                b.ModificationType.Equals("Common Biological") ||
                b.ModificationType.Equals("Metal") ||
                b.ModificationType.Equals("Less Common")
            ).Select(b => (b.ModificationType, b.IdWithMotif)).ToList();
        }

        public List<(string, string)> ListOfModsGptmd { get; set; }
    }
}