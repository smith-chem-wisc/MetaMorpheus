using MassSpectrometry;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
namespace EngineLayer.CrosslinkSearch
{
    public class GlycoPeptides
    {
        public static IEnumerable<Tuple<int, List<Product>>> GlyGetTheoreticalFragments(DissociationType dissociationType,
    List<int> possibleCrosslinkerPositions, PeptideWithSetModifications peptide, Glycan glycan)
        {
            return null;
        }
    }
}
