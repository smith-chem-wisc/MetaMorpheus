using OldInternalLogic;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class ProteinGroup
    {
        public HashSet<Protein> proteins { get; private set; }
        public List<NewPsmWithFDR> psmList { get; private set; }
        public List<CompactPeptide> peptideList { get; private set; }
        public List<CompactPeptide> uniquePeptideList { get; private set; }
        public readonly bool isDecoy;
        private readonly double proteinGroupScore;
        private readonly double summedIntensity;
        private readonly double summedUniquePeptideIntensity;

        public ProteinGroup(HashSet<Protein> proteins, List<NewPsmWithFDR> psmList, List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
        {
            this.proteins = proteins;
            this.psmList = psmList;

            // if any of the proteins in the protein group are decoys, the protein group is a decoy
            foreach (var protein in proteins)
            {
                if (protein.isDecoy)
                    isDecoy = true;
            }

            // build list of compact peptides associated with the protein group
            // all peptides in the group are associated with all proteins in the group
            // if encountering duplicate peptides, only use the best-scoring one
            foreach (var psm in psmList)
            {
                CompactPeptide peptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications);
                peptideList.Add(peptide);
            }

            // calculate the protein group score
            proteinGroupScore = 0;

            // calculate summed psm intensity
            summedIntensity = 0;

            // calculate intensity of only unique peptides
            summedUniquePeptideIntensity = 0;
        }

        public double getScore()
        {
            return proteinGroupScore;
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            // number of proteins in protein group
            foreach (Protein protein in proteins)
                sb.Append("" + protein.FullDescription + " ;; ");
            sb.Append("\t");

            // sequences of proteins in group
            foreach (Protein protein in proteins)
                sb.Append("" + protein.BaseSequence + " ;; ");
            sb.Append("\t");

            // length of each protein
            foreach (Protein protein in proteins)
                sb.Append("" + protein.FullDescription.Length + " ;; ");
            sb.Append("\t");

            // number of proteins in group
            sb.Append("" + proteins.Count);
            sb.Append("\t");

            // number of psm's for the group
            sb.Append("" + psmList.Count);
            sb.Append("\t");

            // number of unique peptides
            sb.Append("" + psmList.Distinct().Count());
            sb.Append("\t");

            // summed psm precursor intensity
            sb.Append(summedIntensity);
            sb.Append("\t");

            // summed unique peptide precursor intensity
            sb.Append(summedUniquePeptideIntensity);
            sb.Append("\t");

            // summed metamorpheus score
            sb.Append(proteinGroupScore);
            sb.Append("\t");

            // isdecoy
            sb.Append(isDecoy);
            sb.Append("\t");

            return sb.ToString();
        }
    }
}