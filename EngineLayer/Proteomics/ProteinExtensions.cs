using Proteomics;
using System.Collections.Generic;

namespace EngineLayer
{
    public static class ProteinExtensions
    {
        /// <summary>
        /// Gets peptides for digestion of a protein
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="digestionParams"></param>
        /// <param name="allKnownFixedModifications"></param>
        /// <param name="variableModifications"></param>
        /// <returns></returns>
        public static IEnumerable<PeptideWithSetModifications> Digest(this Protein protein, IDigestionParams digestionParams, IEnumerable<ModificationWithMass> allKnownFixedModifications, List<ModificationWithMass> variableModifications)
        {
            ProteinDigestion digestion = new ProteinDigestion(digestionParams, allKnownFixedModifications, variableModifications);
            return digestionParams.SemiProteaseDigestion ? digestion.SemiSpecificDigestion(protein) : digestion.Digestion(protein);
        }
    }
}