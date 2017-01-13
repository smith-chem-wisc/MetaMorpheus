using System;
using System.Collections.Generic;
using OldInternalLogic;
using InternalLogicEngineLayer;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InternalLogicEngineLayer
{
    public class ProteinGroup
    {
        public List<Protein> proteins { get; private set; }
        public List<NewPsmWithFDR> psmList { get; private set; }
        public List<CompactPeptide> peptideList { get; private set; }
        public readonly bool isDecoy;
        private double proteinGroupScore = 1;
        private double summedIntensity = 0;
        private double summedUniquePeptideIntensity = 0;

        public ProteinGroup(List<Protein> proteins, List<NewPsmWithFDR> psmList, List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
        {
            this.proteins = proteins;
            this.psmList = psmList;

            // if any of the proteins in the protein group are decoys, the protein group is a decoy
            foreach (var protein in proteins)
            {
                if (protein.isDecoy)
                    this.isDecoy = true;
            }

            // build list of peptides associated with the protein group
            // all peptides in the group are associated with all proteins in the group
            foreach(var psm in psmList)
            {
                peptideList.Add(psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications));
            }

            // find the best scoring psm for each protein in the protein group and multiply them to get the protein group score
            foreach(var protein in proteins)
            {
                double bestScoringPsmScore = 0;

                foreach(var psm in psmList)
                {
                    if (psm.thisPSM.Score > bestScoringPsmScore)
                    {
                        bestScoringPsmScore = psm.thisPSM.Score;
                    }
                }

                proteinGroupScore *= bestScoringPsmScore;
            }
        }

        public double getScore()
        {
            return proteinGroupScore;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            // number of proteins in protein group
            foreach(Protein protein in proteins)
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

            // decoy
            sb.Append(isDecoy);
            sb.Append("\t");

            return sb.ToString();
        }
    }
}
