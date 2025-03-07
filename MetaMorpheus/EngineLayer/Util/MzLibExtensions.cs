using Omics;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer;

public static class MzLibExtensions
{
    // Temporary method unit IBioPolymer has an IsVariant function. 
    public static bool IsVariantPeptide(this IBioPolymerWithSetMods bpwsm)
    {
        if (bpwsm is PeptideWithSetModifications pwsm)
            return pwsm.IsVariantPeptide();
        return false;
    }
}