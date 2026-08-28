using MassSpectrometry;
using Omics.Digestion;
using System.Collections.Generic;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Transcriptomics.Digestion;

namespace EngineLayer;

public static class MzlibExtensions
{
    /// <summary>
    /// Name of the agent that digested the sample, which is what determines the analytes a file can
    /// contain. Use this rather than <see cref="IDigestionParams.DigestionAgent"/> when comparing files:
    /// a fully non-specific search replaces the protease with singleN or singleC, so DigestionAgent
    /// names the search strategy instead of the enzyme, and trypsin and Glu-C files become
    /// indistinguishable. SpecificProtease keeps the enzyme, and the non-specific search preserves it
    /// across the clone it performs. Returns null when unknown.
    /// </summary>
    public static string DigestionAgentName(this IDigestionParams digestionParams)
    {
        if (digestionParams is DigestionParams proteolytic)
        {
            // Unguarded on purpose: DigestionParams sets SpecificProtease unconditionally in its
            // constructor, from a ProteaseDictionary lookup that throws rather than yielding null, and
            // mzLib itself dereferences it bare in ToString and Equals. A null here is therefore not a
            // state this can be reached in, and a fallback for it would be untestable by construction.
            return proteolytic.SpecificProtease.Name;
        }

        return digestionParams?.DigestionAgent?.Name;
    }

    public static Dictionary<DissociationType, List<ProductType>> ProductsFromDissociationType(this IDigestionParams digestionParams)
    {
        if (digestionParams is RnaDigestionParams)
            return Omics.Fragmentation.Oligo.DissociationTypeCollection.ProductsFromDissociationType;
        else
            return Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType;
    }
}
