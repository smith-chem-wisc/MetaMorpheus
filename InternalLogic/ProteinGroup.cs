using System;
using System.Collections.Generic;
using OldInternalLogic;
using InternalLogicEngineLayer;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InternalLogicEngineLayer
{
    class ProteinGroup
    {
        public List<Protein> proteins { get; private set; }
        public List<NewPsmWithFDR> psmList { get; private set; }
        public List<CompactPeptide> peptideList { get; private set; }
        public readonly bool isDecoy;
        private double summedMetaMorpheusScore;
        private double summedIntensity = 0;
        private double summedUniquePeptideIntensity = 0;

        public ProteinGroup(List<Protein> proteins, List<NewPsmWithFDR> psmList, List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
        {
            this.proteins = proteins;
            this.psmList = psmList;
            this.isDecoy = proteins.First().isDecoy;

            foreach(var psm in psmList)
            {
                peptideList.Add(psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications));
                summedMetaMorpheusScore += psm.thisPSM.Score;
            }
        }

        public double getScore()
        {
            return summedMetaMorpheusScore;
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
            sb.Append(summedMetaMorpheusScore);
            sb.Append("\t");

            // decoy
            sb.Append(isDecoy);
            sb.Append("\t");

            return sb.ToString();
        }
    }
}
