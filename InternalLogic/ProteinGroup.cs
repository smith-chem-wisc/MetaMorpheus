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
        public double proteinGroupScore;

        #endregion Public Fields

        #region Internal Constructors

        internal ProteinGroup(HashSet<Protein> proteins, List<NewPsmWithFdr> psmList, HashSet<CompactPeptide> allUniquePeptides, List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
        {
            Proteins = proteins;
            PsmList = psmList;
            PeptideList = new HashSet<CompactPeptide>();
            UniquePeptideList = new HashSet<CompactPeptide>();
            proteinGroupScore = 0;
            QValue = 0;

            // if any of the proteins in the protein group are decoys, the protein group is a decoy
            foreach (var protein in proteins)
            {
                if (protein.IsDecoy)
                    isDecoy = true;
            }

            // build list of compact peptides associated with the protein group
            foreach (var psm in psmList)
            {
                CompactPeptide peptide = psm.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications);
                PeptideList.Add(peptide);
                proteinGroupScore += psm.thisPSM.Score;

                // construct list of unique peptides
                if (allUniquePeptides.Contains(peptide))
                {
                    UniquePeptideList.Add(peptide);
                    //proteinGroupScore += psm.thisPSM.Score;
                }
            }

            if (UniquePeptideList.Count == 0)
                proteinGroupScore = 0;
        }

        internal ProteinGroup(HashSet<Protein> proteins, HashSet<CompactPeptide> peptides, HashSet<CompactPeptide> allUniquePeptides)
        {
            Proteins = proteins;
            PeptideList = peptides;
            UniquePeptideList = new HashSet<CompactPeptide>();
            PsmList = null;
            proteinGroupScore = 0;
            QValue = 0;
            isDecoy = false;

            // if any of the proteins in the protein group are decoys, the protein group is a decoy
            foreach (var protein in proteins)
            {
                if (protein.IsDecoy)
                    isDecoy = true;
            }

            // build list of unique peptides associated with the protein group
            foreach (var peptide in peptides)
            {
                if (allUniquePeptides.Contains(peptide))
                {
                    UniquePeptideList.Add(peptide);
                }
            }
        }

        #endregion Internal Constructors

        #region Public Properties

        public static string TabSeparatedHeader
        {
            get
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
        }

        public HashSet<Protein> Proteins { get; private set; }
        public List<NewPsmWithFdr> PsmList { get; set; }
        public HashSet<CompactPeptide> PeptideList { get; set; }
        public HashSet<CompactPeptide> UniquePeptideList { get; set; }
        public double QValue { get; set; }
        public int cumulativeTarget { get; set; }
        public int cumulativeDecoy { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();

            // list of proteins in the group
            foreach (Protein protein in Proteins)
                sb.Append("" + protein.FullDescription + " ;; ");
            sb.Append("\t");

            // number of proteins in group
            sb.Append("" + Proteins.Count());
            sb.Append("\t");

            // list of unique peptides
            foreach (CompactPeptide peptide in UniquePeptideList)
            {
                string peptideBaseSequence = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                sb.Append("" + peptideBaseSequence + " ;; ");
            }
            sb.Append("\t");

            // list of shared peptides
            foreach (CompactPeptide peptide in PeptideList)
            {
                if (!UniquePeptideList.Contains(peptide))
                {
                    string peptideBaseSequence = string.Join("", peptide.BaseSequence.Select(b => char.ConvertFromUtf32(b)));

                    sb.Append("" + peptideBaseSequence + " ;; ");
                }
            }
            sb.Append("\t");

            // number of unique peptides
            sb.Append("" + UniquePeptideList.Count());
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

        public void scoreThisProteinGroup()
        {
            if (UniquePeptideList.Any())
            {
                foreach (var psm in PsmList)
                {
                    proteinGroupScore += psm.thisPSM.Score;
                }
            }

            else
                proteinGroupScore = 0;
        }

        #endregion Public Methods

    }
}
