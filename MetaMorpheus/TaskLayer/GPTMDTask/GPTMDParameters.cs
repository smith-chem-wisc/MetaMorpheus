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

        /// <summary>
        /// Also try each modification the database put on an identified peptide (a UniProt phosphosite, say) at the
        /// peptide's other residues, and add the sites that score at least as well: on the same amino acid, a tie is
        /// enough; on another, the spectrum must prefer it. The database's own site is kept. Off by default, and
        /// off for any task file that does not name it. On a phospho-enriched set (PXD026608, 12 rat heart files) it
        /// added 1.6% more sites and 0.06% to the database's size.
        /// </summary>
        public bool AddAlternativeSitesForDatabaseMods { get; set; }
    }
}