using Proteomics;

namespace EngineLayer
{
    internal class ProteinInterval
    {
        internal int OneBasedStartResidueInProtein { get; set; }
        internal int OneBasedEndResidueInProtein { get; set; }
        internal Protein Protein { get; set; }
        internal int MissedCleavages { get; set; }
        internal string PeptideString { get; set; }

        internal ProteinInterval(int oneBasedStartResidueInProtein, int oneBasedEndResidueInProtein, Protein protein, int missedCleavages, string peptideString)
        {
            OneBasedStartResidueInProtein = oneBasedStartResidueInProtein;
            OneBasedEndResidueInProtein = oneBasedEndResidueInProtein;
            Protein = protein;
            MissedCleavages = missedCleavages;
            PeptideString = peptideString;
        }
    }
}