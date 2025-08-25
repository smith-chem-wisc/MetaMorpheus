using MassSpectrometry;
using Omics.Digestion;
using System.Collections.Generic;
using Omics.Fragmentation;
using Transcriptomics.Digestion;

namespace EngineLayer;

public static class MzlibExtensions
{
    public static Dictionary<DissociationType, List<ProductType>> ProductsFromDissociationType(this IDigestionParams digestionParams)
    {
        if (digestionParams is RnaDigestionParams)
            return Omics.Fragmentation.Oligo.DissociationTypeCollection.ProductsFromDissociationType;
        else
            return Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType;
    }
}
