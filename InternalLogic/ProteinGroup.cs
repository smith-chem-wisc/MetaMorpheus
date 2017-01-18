using OldInternalLogic;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace InternalLogicEngineLayer
{
    public class ProteinGroup
    {

        #region Public Fields

        public readonly bool isDecoy;
        public readonly double proteinGroupScore;

        #endregion Public Fields

        #region Private Fields

        #endregion Private Fields

        #region Public Constructors

        public ProteinGroup(HashSet<Protein> proteins, List<NewPsmWithFDR> psmList, HashSet<CompactPeptide> allUniquePeptides, List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
        {
            this.proteins = proteins;
            this.psmList = psmList;
            peptideList = new List<CompactPeptide>();
            uniquePeptideList = new List<CompactPeptide>();
            proteinGroupScore = 0;
            QValue = 0;

            // if any of the proteins in the protein group are decoys, the protein group is a decoy
            foreach (var protein in proteins)
            {
                if (protein.isDecoy)
                    isDecoy = true;
            }

            // build list of compact peptides associated with the protein group
            foreach (var psm in psmList)
            {
                CompactPeptide peptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications);
                peptideList.Add(peptide);
                proteinGroupScore += psm.thisPSM.Score;

                // construct list of unique peptides
                if (allUniquePeptides.Contains(peptide))
                {
                    uniquePeptideList.Add(peptide);
                    //proteinGroupScore += psm.thisPSM.Score;
                }
            }

            if (uniquePeptideList.Count == 0)
                proteinGroupScore = 0;
        }

        #endregion Public Constructors

        #region Public Properties

        public readonly HashSet<Protein> proteins;
        public readonly List<NewPsmWithFDR> psmList;
        public readonly List<CompactPeptide> peptideList;
        public readonly List<CompactPeptide> uniquePeptideList;
        public double QValue { get; set; }
        public int cumulativeTarget { get; set; }
        public int cumulativeDecoy { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();

            sb.Append("Proteins in group" + '\t');
            sb.Append("Number of proteins in group" + '\t');
            sb.Append("Unique peptides" + '\t');
            sb.Append("Shared peptides" + '\t');
            sb.Append("Number of unique peptides" + '\t');
            sb.Append("Summed MetaMorpheus Score" + '\t');
            sb.Append("Decoy?" + '\t');
            sb.Append("Cumulative Target" + '\t');
            sb.Append("Cumulative Decoy" + '\t');
            sb.Append("Q-Value (%)");

            return sb.ToString();
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            // list of proteins in the group
            foreach (Protein protein in proteins)
                sb.Append("" + protein.FullDescription + " ;; ");
            sb.Append("\t");

            // number of proteins in group
            sb.Append("" + proteins.Count());
            sb.Append("\t");

            // list of unique peptides
            foreach (CompactPeptide peptide in uniquePeptideList)
            {
                string peptideBaseSequence = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                sb.Append("" + peptideBaseSequence + " ;; ");
            }
            sb.Append("\t");

            // list of shared peptides
            foreach (CompactPeptide peptide in peptideList)
            {
                if (!uniquePeptideList.Contains(peptide))
                {
                    string peptideBaseSequence = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                    sb.Append("" + peptideBaseSequence + " ;; ");
                }
            }
            sb.Append("\t");

            // number of unique peptides
            sb.Append("" + uniquePeptideList.Count());
            sb.Append("\t");

            // summed metamorpheus score
            sb.Append(proteinGroupScore);
            sb.Append("\t");

            // isDecoy
            sb.Append(isDecoy);
            sb.Append("\t");

            // cumulative target
            sb.Append(cumulativeTarget);
            sb.Append("\t");

            // cumulative decoy
            sb.Append(cumulativeDecoy);
            sb.Append("\t");

            // q value
            sb.Append(QValue * 100);
            sb.Append("\t");

            return sb.ToString();
        }

        #endregion Public Methods

    }
}