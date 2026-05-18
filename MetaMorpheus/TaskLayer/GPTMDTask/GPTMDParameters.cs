using System.Collections.Generic;
using System.Linq;
using EngineLayer.Gptmd;

namespace EngineLayer
{
    public class GptmdParameters
    {
        /// <summary>
        /// The default modifications for a G-PTM-D search can be altered here by adding to b.Modification.Type.Equals("  NAME OF MODS TO ADD TO DEFAULT ")
        /// </summary>
        public GptmdParameters()
        {
            GptmdFilters = new();
            ListOfModsGptmd = GlobalVariables.AllModsKnown.Where(b =>
                b.ModificationType.Equals("Common Artifact")
                || b.ModificationType.Equals("Common Biological")
                || b.ModificationType.Equals("Metal") 
                //|| b.ModificationType.Equals("Less Common")
            ).Select(b => (b.ModificationType, b.IdWithMotif)).ToList();
        }

        public List<(string, string)> ListOfModsGptmd { get; set; }
        public List<IGptmdFilter> GptmdFilters { get; set; } 
        public bool WriteDecoys { get; set; } = true;
    }
}