using Proteomics;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class ProteinDigestion
    {
        public Protease Protease { get; set; }
        public int MaximumMissedCleavages { get; set; }
        public IDigestionParams DigestionParams { get; set; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }
        public int? MinPeptidesLength { get; set; }
        public int? MaxPeptidesLength { get; set; }
        public IEnumerable<ModificationWithMass> AllKnownFixedModifications { get; set; }
        public List<ModificationWithMass> VariableModifications { get; set; }

        /// <summary>
        /// Initializes digestion object
        /// </summary>
        /// <param name="digestionParams"></param>
        /// <param name="allKnownFixedModifications"></param>
        /// <param name="variableModifications"></param>
        public ProteinDigestion(IDigestionParams digestionParams, IEnumerable<ModificationWithMass> allKnownFixedModifications, List<ModificationWithMass> variableModifications)
        {
            DigestionParams = digestionParams;
            Protease = digestionParams.Protease;
            MaximumMissedCleavages = digestionParams.MaxMissedCleavages;
            InitiatorMethionineBehavior = digestionParams.InitiatorMethionineBehavior;
            MinPeptidesLength = digestionParams.MinPeptideLength;
            MaxPeptidesLength = digestionParams.MaxPeptideLength;
            AllKnownFixedModifications = allKnownFixedModifications;
            VariableModifications = variableModifications;
        }

        /// <summary>
        /// Gets peptides for semispecific digestion of a protein
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        public IEnumerable<PeptideWithSetModifications> SemiSpecificDigestion(Protein protein)
        {
            List<Peptide> intervals = new List<Peptide>();
            List<int> oneBasedIndicesToCleaveAfter = Protease.GetDigestionSiteIndices(protein.BaseSequence);

            for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - MaximumMissedCleavages - 1; i++)
            {
                if (Protease.Retain(i, InitiatorMethionineBehavior, protein[0])
                    && Protease.OkayLength(oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - oneBasedIndicesToCleaveAfter[i], MinPeptidesLength, MaxPeptidesLength))
                {
                    intervals.Add(new Peptide(protein, oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1],
                        oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - oneBasedIndicesToCleaveAfter[i], "semi"));
                }

                if (Protease.Cleave(i, InitiatorMethionineBehavior, protein[0])
                    && Protease.OkayLength(oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - 1, MinPeptidesLength, MaxPeptidesLength))
                {
                    intervals.Add(new Peptide(protein, 2, oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1],
                        oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - 1, "semi:M cleaved"));
                }
            }

            int lastIndex = oneBasedIndicesToCleaveAfter.Count - 1;
            int maxIndex = MaximumMissedCleavages < lastIndex ? MaximumMissedCleavages : lastIndex;
            for (int i = 1; i <= maxIndex; i++)
            {
                if (DigestionParams.TerminusTypeSemiProtease == TerminusType.N) //tricky, it's N because we want the extra peptide at the C terminus |_
                {
                    if (Protease.OkayLength(oneBasedIndicesToCleaveAfter[lastIndex] - oneBasedIndicesToCleaveAfter[lastIndex - i], MinPeptidesLength, MaxPeptidesLength))
                    {
                        intervals.Add(new Peptide(protein, oneBasedIndicesToCleaveAfter[lastIndex - i] + 1, oneBasedIndicesToCleaveAfter[lastIndex],
                            oneBasedIndicesToCleaveAfter[lastIndex] - oneBasedIndicesToCleaveAfter[lastIndex - i], "semiN"));
                    }
                }
                else //TerminusType.C
                {
                    if (Protease.OkayLength(oneBasedIndicesToCleaveAfter[i] - oneBasedIndicesToCleaveAfter[0], MinPeptidesLength, MaxPeptidesLength))
                    {
                        intervals.Add(new Peptide(protein, oneBasedIndicesToCleaveAfter[0] + 1, oneBasedIndicesToCleaveAfter[i],
                            oneBasedIndicesToCleaveAfter[i] - oneBasedIndicesToCleaveAfter[0], "semiC"));
                    }
                }
            }

            // Also digest using the proteolysis product start/end indices
            intervals.AddRange(
                protein.ProteolysisProducts
                    .Where(proteolysisProduct => proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length)
                    .Select(proteolysisProduct => new Peptide(protein, proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value,
                        0, proteolysisProduct.Type + " start")));

            return intervals.SelectMany(peptide => peptide.GetModifiedPeptides(AllKnownFixedModifications, DigestionParams, VariableModifications));
        }

        /// <summary>
        /// Gets peptides for specific protease digestion of a protein
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        public IEnumerable<PeptideWithSetModifications> Digestion(Protein protein)
        {
            var intervals = Protease.GetDigestionIntervals(protein, MaximumMissedCleavages, InitiatorMethionineBehavior, MinPeptidesLength, MaxPeptidesLength);
            return intervals.SelectMany(peptide => peptide.GetModifiedPeptides(AllKnownFixedModifications, DigestionParams, VariableModifications));
        }
    }
}